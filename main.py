# -*- coding: utf-8 -*-
"""
Created on Thu May 23 14:30:56 2019

@author: ryoichi.tsurumaki

このスクリプトは、集中型流出解析モデルを用いて河川水位予測を
行うためのものです。
"""

import pprint
import subprocess
import os
import sys
import spotpy
#from spotpy import examples
#from spotpy.examples importp spot_setup_rosenbrock
from spotpy.examples.spot_setup_rosenbrock import spot_setup
from spotpy import analyser
import datetime
from tankModelCalculator import TankModelCalculator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from paramOptimizer import ParamOptimizer
from kalmanCalculator import KalmanCalculator
from scipy.optimize import differential_evolution
from dataReader import DataReader
from graphMaker import GraphMaker


if __name__ == '__main__':
    print(__doc__)
    startCalTime = datetime.datetime.now()

    # データ入力
    # ----------------------------------------
    allDataDF = pd.DataFrame()
    dataReader = DataReader(allDataDF)
    cal_settings = dataReader.getCalSettings()
    dataReader.putTime(cal_settings["startTime"],
                       cal_settings["endTime"],
                       cal_settings["timeInterval"],
                       cal_settings["timescale"])
    dataReader.setInputFilePath()
    dataReader.readRainfallData(cal_settings["used_rainfallData"],
                                cal_settings["timeInterval"],
                                cal_settings["timescale"],
                                cal_settings["startTime"],
                                cal_settings["endTime"])
    dataReader.readWaterLevelData(cal_settings["used_waterLevelData"],
                                  cal_settings["timeInterval"],
                                  cal_settings["timescale"],
                                  cal_settings["startTime"],
                                  cal_settings["endTime"])
    HQParams = dataReader.getHQParams()
    allDataDF = dataReader.getInputDataDF()

    # HQ 式から流量の計算
    # ----------------------------------------
    for x in HQParams:
        allDataDF.loc[(allDataDF["water level(obs)"] >= x["H_min"]) &
                      (allDataDF["water level(obs)"] < x["H_max"]),
                      "flow rate(HQ)"] = \
                      x["a"] * (allDataDF["water level(obs)"] - x["b"])**2
    allDataDF.loc[allDataDF["water level(obs)"] < HQParams[0]["H_min"],
                  "flow rate(HQ)"] = 0  # for less than the minimum water level
    allDataDF.loc[allDataDF["water level(obs)"] >= HQParams[-1]["H_max"],
                  "flow rate(HQ)"] = np.nan

    # 各種計算
    # ----------------------------------------
    totalRainfall = allDataDF["rainfall"].loc[
        cal_settings["rainfallStartTime"] + datetime.timedelta(hours=1):
        cal_settings["rainfallEndTime"]].sum()  # 総雨量
    totalRainfallExceptLoss = allDataDF["rainfall"].loc[
        cal_settings["floodStartTime"] + datetime.timedelta(hours=1):
        cal_settings["rainfallEndTime"]].sum()  # 初期損失雨量をのぞく総雨量
    # 初期損失雨量
    initialLossRainfall = totalRainfall - totalRainfallExceptLoss
    # 基底流出高(一定を仮定=0)
    baseFlowHeight = allDataDF.at[cal_settings["floodStartTime"], "flow rate(HQ)"]
    # 直接流出高のデータフレーム
    directFlowHeightDF = allDataDF["flow rate(HQ)"].loc[
        cal_settings["floodStartTime"]:cal_settings["floodEndTime"]] - baseFlowHeight
    # 総直接流出高
    totalDirectFlowHeight = directFlowHeightDF.sum()
    # 流出率
    flowRatio = totalDirectFlowHeight / totalRainfallExceptLoss
    # 有効雨量のデータフレーム
    allDataDF["effective rainfall"] = list(allDataDF["rainfall"] * flowRatio)
    # effectiveRainfallDF = allDataDF["rainfall"].\
    # loc[flowStartTime+datetime.timedelta(hours=1):flowEndTime] * flowRatio
    # effectiveRainfallDF = allDataDF["rainfall"] * flowRatio

    # 流出計算
    # ----------------------------------------
    if cal_settings["used_algorithm"] == "kalmanFilter":
        kalmanCalculator = KalmanCalculator(cal_settings)
        kalman_settings = dataReader.getKalmanSettings()
        kalmanCalculator.setReadData(allDataDF, cal_settings["catchmentArea"],
                                     kalman_settings)
        kalmanCalculator.run_EKF()
    elif cal_settings["used_algorithm"] == "de":
        paramOptimizer = ParamOptimizer(cal_settings)
        tankModel_settings = dataReader.getTankParams()
        paramOptimizer.setReadDataForTank(allDataDF,
                                          cal_settings["catchmentArea"],
                                          tankModel_settings)
        bounds = paramOptimizer.getBounds()
        result = differential_evolution(paramOptimizer.searchingFunc, bounds)
        # 最適化されたパラメータを用いた流出モデルの計算
        # modelError = paramOptimizer.searchingFunc(result["x"])

    # 結果の整理
    # ----------------------------------------
    def summarizeResults(headerList):
        outputDF = pd.DataFrame(outputList, columns=headerList)
        outputDF = outputDF.set_index("Date")
        outputDF.loc[outputDF["flow rate"] < HQParams[0]["Q_min"],
                     "water level"] = np.nan
        # 水位計算(HQ 式の適用範囲に注意)
        # ----------------------------------------
        outputDF.loc[(outputDF["flow rate"] >= HQParams[0]["Q_min"]) & \
                     (outputDF["flow rate"] < HQParams[0]["Q_max"]), "water level"] =\
                      HQParams[0]["b"] + (outputDF["flow rate"]/HQParams[0]["a"])**0.5
        outputDF.loc[(outputDF["flow rate"] >= HQParams[0]["Q_max"]) & \
                     (outputDF["flow rate"] < HQParams[1]["Q_max"]), "water level"] =\
                      HQParams[1]["b"] + (outputDF["flow rate"]/HQParams[1]["a"])**0.5
        outputDF.loc[(outputDF["flow rate"] >= HQParams[1]["Q_max"]) &\
                     (outputDF["flow rate"] < HQParams[2]["Q_max"]), "water level"] =\
                      HQParams[2]["b"] + (outputDF["flow rate"]/HQParams[2]["a"])**0.5
        outputDF.loc[outputDF["flow rate"] > HQParams[2]["Q_max"], "water level"] = np.nan
        allDataDF["flow rate(cal)"] = outputDF["flow rate"]
        allDataDF["water level(cal)"] = outputDF["water level"]
        return outputDF

    outputFilePath = dataReader.getOutputFilePath()
    if cal_settings["used_algorithm"] == "de":
        outputList = paramOptimizer.getOutputData()
        if cal_settings["used_flowModel"] == "tankModel":
            headerList = ["storage height", "tn1", "tn2", "tn3", "total flow height",
                          "flow height1", "flow height2", "flow height3",
                          "leakage1", "leakage2", "leakage3", "Date",
                          "flow rate", "rainfall"]
            outputDF = summarizeResults(headerList)
            allDataDF["storage height"] = outputDF["storage height"]
        elif cal_settings["used_flowModel"] == "classicOneValueStorageFunc" or \
            cal_settings["used_flowModel"] == "classicTwoValueStorageFunc":
            headerList = ["Date", "flow rate", "rainfall", "storage height"]
            outputDF = summarizeResults(headerList)
            allDataDF["storage height"] = outputDF["storage height"]
        elif cal_settings["used_flowModel"] == "twoStepTwoValueStorageFunc":
            headerList = ["Date", "flow rate", "rainfall", "storage height"]
            outputDF = summarizeResults(headerList)
            allDataDF["storage height"] = outputDF["storage height"]
    elif cal_settings["used_algorithm"] == "kalmanFilter":
        outputList = kalmanCalculator.getOutputData()
        headerList = ["Date", "flow rate", "rainfall", "flow rate error"]
        outputDF = summarizeResults(headerList)
        # allDataDF["storage height"] = outputDF["storage height"]

    # 図の出力
    # ----------------------------------------
    graphMaker = graphMaker(allDataDF)

    elapsedTime = (datetime.datetime.now() - startCalTime).total_seconds()
    print(f"##################")
    print(f"** CALCULATE TIME: {int(elapsedTime//60)}min {int(elapsedTime%60)}sec")
    print(f"##################")
