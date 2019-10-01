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
#from spotpy.examples.spot_setup_rosenbrock import spot_setup
#from spotpy import analyser
import datetime
from tankModelCalculator import TankModelCalculator
import pandas as pd
import numpy as np
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

    if cal_settings["used_rainfallData"]:
        # 雨量 (mm/hr)
        dataReader.readRainfallData(cal_settings["used_rainfallData"],
                                    cal_settings["timeInterval"],
                                    cal_settings["timescale"],
                                    cal_settings["startTime"],
                                    cal_settings["endTime"])
    if cal_settings["used_waterLevelData"]:
        # 水位 (m)
        dataReader.readWaterLevelData(cal_settings["used_waterLevelData"],
                                      cal_settings["timeInterval"],
                                      cal_settings["timescale"],
                                      cal_settings["startTime"],
                                      cal_settings["endTime"])
        # HQ 式
        HQParams = dataReader.getHQParams()
    if cal_settings["used_flowRateData"]:
        # 流量 (m^3/s)
        dataReader.readFlowRateData(cal_settings["used_flowRateData"],
                                    cal_settings["timeInterval"],
                                    cal_settings["timescale"],
                                    cal_settings["startTime"],
                                    cal_settings["endTime"])
    # 流域面積 (km^2)
    if cal_settings["obsName"] == "atsumaohashi":
        catchmentArea = 238.4
    elif cal_settings["obsName"] == "motsukisamu":
        catchmentArea = 6.13

    allDataDF = dataReader.getInputDataDF()

    # HQ 式から流量の計算
    # ----------------------------------------
    if cal_settings["used_waterLevelData"]:
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
    baseFlowRate = allDataDF["flow rate(HQ)"].min()  # 基底流出高(期間平均)
    totalRainfall = allDataDF["rainfall"].sum()  # 総雨量
    # 直接流出高
    directFlowVolumeDF = (allDataDF["flow rate(HQ)"] - baseFlowRate) * \
        (cal_settings["endTime"] - cal_settings["startTime"]).total_seconds()
    totalDirectFlowVolume = directFlowVolumeDF.sum()  # 総直接流出高
    # 流出率
    flowRatio = totalDirectFlowVolume / \
                (totalRainfall/1000 * catchmentArea * 10**6)
    if flowRatio > 1:
        flowRatio = 1
    # 有効雨量のデータフレーム
    allDataDF["effective rainfall"] = list(allDataDF["rainfall"] * flowRatio)
    """
    totalRainfallExceptLoss = allDataDF["rainfall"].loc[
        cal_settings["floodStartTime"] + datetime.timedelta(hours=1):
        cal_settings["rainfallEndTime"]].sum()  # 初期損失雨量をのぞく総雨量
    # 初期損失雨量
    initialLossRainfall = totalRainfall - totalRainfallExceptLoss
    # 直接流出高のデータフレーム
    directFlowHeightDF = allDataDF["flow rate(HQ)"].loc[
        cal_settings["floodStartTime"]:cal_settings["floodEndTime"]] - baseFlowHeight
    # 総直接流出高
    totalDirectFlowHeight = directFlowHeightDF.sum()
    # 流出率
    flowRatio = totalDirectFlowHeight / totalRainfallExceptLoss
    # effectiveRainfallDF = allDataDF["rainfall"].\
    # loc[flowStartTime+datetime.timedelta(hours=1):flowEndTime] * flowRatio
    # effectiveRainfallDF = allDataDF["rainfall"] * flowRatio
    """

    # 流出計算
    # ----------------------------------------
    if cal_settings["used_algorithm"] == "kalmanFilter":
        kalmanCalculator = KalmanCalculator(cal_settings)
        kalman_settings = dataReader.getKalmanSettings()
        kalmanCalculator.setReadOnlyData(allDataDF, catchmentArea, kalman_settings,
                                         baseFlowRate, flowRatio)
        kalmanCalculator.run_EKF()
    elif cal_settings["used_algorithm"] == "de":
        paramOptimizer = ParamOptimizer(cal_settings, baseFlowRate, catchmentArea, allDataDF)
        if cal_settings["used_flowModel"] == "tankModel":
            tankModel_settings = dataReader.getTankParams()
            paramOptimizer.setReadDataForTank(tankModel_settings)
        bounds = paramOptimizer.getBounds()
        resultTuple = differential_evolution(paramOptimizer.searchFunc, bounds)
        # result = resultTuple[0]
        #populationList = resultTuple[1]
        #population_energyList = resultTuple[2]
        # 最適化されたパラメータを用いた流出モデルの計算
        # modelError = paramOptimizer.searchFunc(result["x"])

    # アニメーションの作成
    # ----------------------------------------
    # graphMaker = GraphMaker()
    # graphMaker.makeAnimation(populationList, population_energyList)

    # 結果の整理
    # ----------------------------------------
    def summarizeResults(headerList, allDataDF, outputList, cal_settings):
        outputDF = pd.DataFrame(outputList, columns=headerList)
        outputDF = outputDF.set_index("Date")
        if cal_settings["used_waterLevelData"]:
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
            allDataDF["water level(cal)"] = outputDF["water level"]
        allDataDF = pd.concat([allDataDF, outputDF["storage height"]], axis=1)
        allDataDF["flow rate(cal)"] = outputDF["flow rate"]
        if cal_settings["used_algorithm"] == "kalmanFilter":
            allDataDF["flow rate error"] = outputDF["flow rate error"]
            allDataDF["rainfall"] = outputDF["rainfall"]
        return allDataDF

    outputFilePath = dataReader.getOutputFilePath()
    if cal_settings["used_algorithm"] == "de":
        outputList = paramOptimizer.getOutputData()
        if cal_settings["used_flowModel"] == "tankModel":
            headerList = ["storage height", "tn1", "tn2", "tn3", "total flow height",
                          "flow height1", "flow height2", "flow height3",
                          "leakage1", "leakage2", "leakage3", "Date",
                          "flow rate", "rainfall"]
            allDataDF = summarizeResults(headerList, allDataDF, outputList, cal_settings)
        elif cal_settings["used_flowModel"] == "classicOneValueStorageFunc" or \
            cal_settings["used_flowModel"] == "classicTwoValueStorageFunc":
            headerList = ["Date", "flow rate", "rainfall", "storage height"]
            allDataDF = summarizeResults(headerList, allDataDF, outputList, cal_settings)
        elif cal_settings["used_flowModel"] == "twoStepTwoValueStorageFunc":
            headerList = ["Date", "flow rate", "rainfall", "storage height"]
            allDataDF = summarizeResults(headerList, allDataDF, outputList, cal_settings)
    elif cal_settings["used_algorithm"] == "kalmanFilter":
        outputList = kalmanCalculator.getOutputData()
        headerList = ["Date", "flow rate", "rainfall", "flow rate error",
                      "storage height"]
        allDataDF = summarizeResults(headerList, allDataDF, outputList, cal_settings)

    # 図の出力
    # ----------------------------------------
    graphMaker = GraphMaker(allDataDF)
    endTime = cal_settings["endTime"] + \
              datetime.timedelta(hours=int(cal_settings["forecastTime"]))
    if cal_settings["used_algorithm"] == "de":
        if cal_settings["used_flowModel"] == "classicOneValueStorageFunc":
            savePath = outputFilePath + \
                       "differential_evolution/oneValueStorageFunction/"
        elif cal_settings["used_flowModel"] == "classicTwoValueStorageFunc":
            savePath = outputFilePath + \
                       "differential_evolution/twoValueStorageFunction/"
        elif cal_settings["used_flowModel"] == "twoStepTwoValueStorageFunc":
            savePath = outputFilePath + \
                       "differential_evolution/twoStepTankStorageFunction/"
        elif cal_settings["used_flowModel"] == "tankModel":
            savePath = outputFilePath + "differential_evolution/tankModel/"
    elif cal_settings["used_algorithm"] == "kalmanFilter":
        savePath = outputFilePath + "kalmanFilter/"
    xTickList = [cal_settings["startTime"], endTime, cal_settings["timescale"]]
    xLabel = "Date"

    # flow rate - date
    yNameList = ["flow rate(HQ)", "flow rate(cal)", "rainfall"]
    if cal_settings["timescale"] == "hours" or "minutes":
        yLabelList = ["Flow rate (m$^3$/s)", "Rainfall (mm/hr)"]
    if cal_settings["timescale"] == "days":
        yLabelList = ["Flow rate (m$^3$/s)", "Rainfall (mm/day)"]
    outputFileName = cal_settings["obsName"] + "_" + \
        cal_settings["used_flowModel"] + "_" + \
        cal_settings["used_algorithm"] + "_" + \
        cal_settings["startTime"].strftime("%Y%m%d%H%M") + "_" + \
        cal_settings["endTime"].strftime("%Y%m%d%H%M") + "_" + \
        cal_settings["timeInterval"] + \
        cal_settings["timescale"] + "_" + "FR.png"
    graphMaker.makeFlowRateGraph(yNameList, xTickList, xLabel, yLabelList)
    # graphMaker.save(outputFileName, savePath)

    # water level - date
    if cal_settings["used_waterLevelData"]:
        yNameList = ["water level(obs)", "water level(cal)", "rainfall"]
        yLabelList = ["Water lebel (m)", "Rainfall (mm/hr)"]
        outputFileName = cal_settings["obsName"] + "_" + \
            cal_settings["used_flowModel"] + "_" + \
            cal_settings["used_algorithm"] + "_" + \
            cal_settings["startTime"].strftime("%Y%m%d%H%M") + "_" + \
            cal_settings["endTime"].strftime("%Y%m%d%H%M") + "_" + \
            cal_settings["timeInterval"] + \
            cal_settings["timescale"] + "_" + "WL.png"
        graphMaker.makeWaterLevelGraph(yNameList, xTickList, xLabel, yLabelList)
        # graphMaker.save(outputFileName, savePath)

    # flow rate - storage height
    xName = "flow rate(cal)"
    yName = "storage height"
    xLabel = "Flow rate (m$^3$/s)"
    yLabel = "Storage height(mm)"
    outputFileName = cal_settings["obsName"] + "_" + \
        cal_settings["used_flowModel"] + "_" + \
        cal_settings["used_algorithm"] + "_" + \
        cal_settings["startTime"].strftime("%Y%m%d%H%M") + "_" + \
        cal_settings["endTime"].strftime("%Y%m%d%H%M") + "_" + \
        cal_settings["timeInterval"] + \
        cal_settings["timescale"] + "_" + "FS.png"
    graphMaker.makeFlow_StorageRelation(xName, yName, xLabel, yLabel)
    # graphMaker.save(outputFileName, savePath)

    elapsedTime = (datetime.datetime.now() - startCalTime).total_seconds()
    print(f"################################")
    print(f"   CALCULATE TIME: {int(elapsedTime//60)}min {int(elapsedTime%60)}sec")
    print(f"################################")
