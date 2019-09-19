# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:24:00 2019

@author: ryoichi.tsurumaki

"""

import numpy as np
import datetime
import pandas as pd
import pprint
from grib2Reader import Grib2Reader
import configparser


class DataReader(object):
    """
    class for read and storage of meteorological data
    """
    def __init__(self, allDataDF=None):
        """
        ----
        Input
        ----
        allDataDF: DataFrame
            * storage place for all data
        inifile: configparser
            * file storaged for the setting
        """
        if allDataDF is not None:
            self.allDataDF = allDataDF
        else:
            self.allDataDF = pd.DataFrame()
        self.inifile = configparser.ConfigParser()
        self.inifile.read('./config.ini', 'utf-8')

    def getInputDataDF(self):
        return self.allDataDF

    def getGrib2Params(self):
        grib2_settings = \
            {"timescale": self.inifile.get("grib2TimeInterval", "timescale"),
             "timeInterval": self.inifile.get("grib2TimeInterval", "timeInterval"),
             "num": self.inifile.get("grib2TimeInterval", "num"),
             "latitude": self.inifile.get("center_location", "latitude"),
             "longitude": self.inifile.get("center_location", "longitude"),
             }

        return grib2_settings

    def getCalSettings(self):
        """
        method for getting fundamental parameters for the calculation
        ----
        output
        ----
        cal_settings: dic
        """
        cal_settings = \
            {"startTime": datetime.datetime.strptime(self.inifile.get(
                "inputFile", "startTime"), "%Y/%m/%d %H:%M:%S"),
             "endTime": datetime.datetime.strptime(self.inifile.get(
                "inputFile", "endTime"), "%Y/%m/%d %H:%M:%S"),
             "timeInterval": self.inifile.get("inputFile", "timeInterval"),
             "timescale": self.inifile.get("inputFile", "timescale"),
             "used_rainfallData": self.inifile.get("inputFile",
                                                   "used_rainfallData"),
             "used_waterLevelData": self.inifile.get("inputFile",
                                                     "used_waterLevelData"),
             "used_algorithm": self.inifile.get("optimizationParameters",
                                                "algorithm"),
             "used_objFunc": self.inifile.get("optimizationParameters",
                                              "objFunction"),
             "used_flowModel": self.inifile.get("optimizationParameters",
                                                "flowModel"),
             "catchmentArea": float(self.inifile.get("riverParameters",
                                                     "catchmentArea")),
             "rainfallStartTime": datetime.datetime.strptime(
                 self.inifile.get("riverParameters", "rainfallStartTime"),
                 "%Y/%m/%d %H:%M:%S"),  # 降雨開始時刻
             "rainfallEndTime": datetime.datetime.strptime(
                 self.inifile.get("riverParameters", "rainfallEndTime"),
                 "%Y/%m/%d %H:%M:%S"),  # 降雨終了時刻
             "floodStartTime": datetime.datetime.strptime(
                 self.inifile.get("riverParameters", "floodStartTime"),
                 "%Y/%m/%d %H:%M:%S"),  # 氾濫開始時刻
             "floodEndTime": datetime.datetime.strptime(
                 self.inifile.get("riverParameters", "floodEndTime"),
                 "%Y/%m/%d %H:%M:%S"),  # 氾濫終了時刻
             "obsDateTime": datetime.datetime.strptime(
                 self.inifile.get("riverParameters", "obsDateTime"),
                 "%Y/%m/%d %H:%M:%S"),  # 観測時刻
             "riverName": self.inifile.get("riverParameters", "riverName"),
             "forecastTime": self.inifile.get("riverParameters", "forecastTime"),
             }

        return cal_settings

    def putTime(self, startTime, endTime, timeInterval, timescale):
        """
        method for putting data about time to the dataframe
        ----
        Input
        ----
        """
        timestamp = startTime
        dateList = [startTime]
        if timescale == "hours":
            dataNum = int((endTime - startTime).
                          total_seconds()/3600/int(timeInterval) + 1)
            gridTimeList = [datetime.timedelta(hours=1)]
        elif timescale == "minutes":
            dataNum = int((endTime - startTime).
                          total_seconds()/60/int(timeInterval) + 1)
            gridTimeList = [datetime.timedelta(minutes=10)]
        # for i in range(dataNum-1):
        for i in range(dataNum-1):
            timestamp += eval("datetime.timedelta(" + timescale + "=" + timeInterval + ")")
            dateList.append(timestamp)
            gridTimeList.append(timestamp - dateList[i])
        # gridTimeList.append(np.nan)
        self.allDataDF["grid time"] = gridTimeList
        self.allDataDF.index = dateList
        # pList = []
        # for delta in gridTimeList:
            # tmpList.append(delta.total_seconds()/3600)
        # tmpList.append(np.nan)
        # self.allDataDF["grid time"] = tmpList
        # self.allDataDF.index = tmpDateList

    def setInputFilePath(self, *args):
        if args:
            self.inputFilePath = args[0]
        else:
            self.inputFilePath = self.inifile.get("inputFile", "path")

    def getOutputFilePath(self, *args):
        if args:
            return args[0]
        else:
            return self.inifile.get("outputFile", "path")

    def getKalmanSettings(self):
        """
        method for getting fundamental parameters for the kalman filter
        ----
        output
        ----
        kalman_settings: dic
        """
        kalman_settings =\
            {"dim_x": int(self.inifile.get("kalmanParameters",
                                           "dimensionOfState")),
             "dim_z": int(self.inifile.get("kalmanParameters",
                                           "dimensionOfObs")),
             "dim_u": int(self.inifile.get("kalmanParameters",
                                           "dimensionOfControl")),
             "p1": float(self.inifile.get("kalmanParameters",
                                          "storageIndex1")),
             "p2": float(self.inifile.get("kalmanParameters",
                                          "storageIndex2")),
             }
        return kalman_settings

    def getTankParams(self):
        """
        method for getting data about the tank model parameters recommended
        by Meteorological Agency
        ----
        Input
        ----
        """
        return \
            {
            # 第 1 タンクの下部横孔の流出係数(hr^-1)
            "alpha1": float(self.inifile.get("tankParameters", "alpha1")),
            # 第 1 タンクの上部横孔の流出係数(hr^-1)
            "alpha2": float(self.inifile.get("tankParameters", "alpha2")),
            # 第 2 タンクの横孔の流出係数(hr^-1)
            "alpha3": float(self.inifile.get("tankParameters", "alpha3")),
            # 第 3 タンクの横孔の流出係数(hr^-1)
            "alpha4": float(self.inifile.get("tankParameters", "alpha4")),
            # 第 1 タンクの下孔の浸透係数(hr^-1)
            "beta1": float(self.inifile.get("tankParameters", "beta1")),
            # 第 1 タンクの下孔の浸透係数(hr^-1)
            "beta2": float(self.inifile.get("tankParameters", "beta2")),
            # 第 2 タンクの下孔の浸透係数(hr^-1)
            "beta3": float(self.inifile.get("tankParameters", "beta3")),
            # 第 1 タンクの下部横孔の高さ(mm)
            "height1": float(self.inifile.get("tankParameters", "height1")),
            # 第 1 タンクの上部横孔の高さ(mm)
            "height2": float(self.inifile.get("tankParameters", "height2")),
            # 第 2 タンクの横孔の高さ(mm)
            "height3": float(self.inifile.get("tankParameters", "height3")),
            # 第 3 タンクの横孔の高さ(mm)
            "height4": float(self.inifile.get("tankParameters", "height4")),
            # 第 1 タンクの初期貯留高(mm)
            "storage1": float(self.inifile.get("tankParameters", "storage1")),
            # 第 2 タンクの初期貯留高(mm)
            "storage2": float(self.inifile.get("tankParameters", "storage2")),
            # 第 3 タンクの初期貯留高(mm)
            "storage3": float(self.inifile.get("tankParameters", "storage3")),
            }

    def getHQParams(self):
        """
        method for getting parameters for the HQ equation
        ----
        """
        H_min1 = float(self.inifile.get("HQparameters", "H_min1"))
        H_max1 = float(self.inifile.get("HQparameters", "H_max1"))
        a1 = float(self.inifile.get("HQparameters", "a1"))
        b1 = float(self.inifile.get("HQparameters", "b1"))
        Q_min1 = a1 * (H_min1 - b1)**2
        Q_max1 = a1 * (H_max1 - b1)**2
        H_min2 = float(self.inifile.get("HQparameters", "H_min2"))
        H_max2 = float(self.inifile.get("HQparameters", "H_max2"))
        a2 = float(self.inifile.get("HQparameters", "a2"))
        b2 = float(self.inifile.get("HQparameters", "b2"))
        Q_min2 = a2 * (H_min2 - b2)**2
        Q_max2 = a2 * (H_max2 - b2)**2
        H_min3 = float(self.inifile.get("HQparameters", "H_min3"))
        H_max3 = float(self.inifile.get("HQparameters", "H_max3"))
        a3 = float(self.inifile.get("HQparameters", "a3"))
        b3 = float(self.inifile.get("HQparameters", "b3"))
        Q_min3 = a3 * (H_min3 - b3)**2
        Q_max3 = a3 * (H_max3 - b3)**2
        range1 =\
            {
            "H_min": H_min1,
            "H_max": H_max1,
            "a": a1,
            "b": b1,
            "Q_min": Q_min1,
            "Q_max": Q_max1,
            }
        range2 =\
            {
            "H_min": H_min2,
            "H_max": H_max2,
            "a": a2,
            "b": b2,
            "Q_min": Q_min2,
            "Q_max": Q_max2,
            }
        range3 =\
            {
            "H_min": H_min3,
            "H_max": H_max3,
            "a": a3,
            "b": b3,
            "Q_min": Q_min3,
            "Q_max": Q_max3,
            }
        return [range1, range2, range3]

    def readRainfallData(self, used_rainfallData, timeInterval, timescale,
                         startTime, endTime):
        """
        method for setting data about the rainfall at the motsukisamu river
        ----
        Input
        ----
        """
        # 札幌アメダスの気象庁データ (1 時間)
        if used_rainfallData == "sapporo amedas" and \
           timeInterval+timescale == "1hours":
            inputFileName = self.inifile.get("inputFile", "rainfallFileName1")
            inputDF = pd.read_csv(self.inputFilePath + inputFileName,
                encoding=None, skiprows=5, header=None, sep=",",
                skipinitialspace=True, usecols=[0, 1], index_col=0)
            inputDF.index = pd.to_datetime(inputDF.index)
            inputDF.columns = ["rainfall"]
            self.allDataDF = pd.concat([self.allDataDF, inputDF], axis=1)

        # 望月寒川の web データ (1 時間)
        elif used_rainfallData == "motsukisamu" and \
             timeInterval+timescale == "1hours":
            # データが格納された列番号(時間ごと)
            dataColumnNo = [x for x in range(7, 31)]
            inputFileName = self.inifile.get("inputFile", "rainfallFileName2")
            inputDF = pd.read_csv(self.inputFilePath + inputFileName,
                                  encoding="shift-jis",
                                  skiprows=1,
                                  header=None,
                                  sep=",",
                                  skipinitialspace=True,
                                  usecols=dataColumnNo,
                                  index_col=None)
            # convert into NaN except for the numerical value
            # inputDF = inputDF.convert_objects(convert_numeric=True)
            inputDF = inputDF.apply(pd.to_numeric, args=('coerce',))
            startYear = int(startTime.strftime("%Y"))
            timestamp = datetime.datetime(startYear, 1, 1, 1)  # from YYYY/1/1 1:00
            dateList = []
            readingRowNo = 0
            rainfallList = []
            while timestamp < datetime.datetime(startYear+1, 1, 1, 0):
                rainfallList.append(list(inputDF.iloc[readingRowNo, 0:24]))
                for i in range(24):
                    dateList.append(timestamp)
                    timestamp += eval("datetime.timedelta(hours=1)")
                readingRowNo += 1
            rainfallList = np.array(rainfallList).flatten()
            tmpSeries = pd.Series(rainfallList, index=dateList)
            self.allDataDF["rainfall"] = tmpSeries[startTime:endTime]

        # 望月寒川の A ロガーデータ (10分)
        elif used_rainfallData == "motsukisamu" and \
           timeInterval+timescale == "10minutes":
            inputFileName = self.inifile.get("inputFile", "rainfallFileName3")
            inputDF = pd.read_csv(self.inputFilePath + inputFileName,
                encoding=None, skiprows=None, header=None, sep=",",
                skipinitialspace=True, usecols=[0, 2], index_col=0)
            inputDF.index = pd.to_datetime(inputDF.index)
            inputDF.columns = ["rainfall"]
            self.allDataDF = pd.concat([self.allDataDF, inputDF/10*6], axis=1,
                                       join_axes=[self.allDataDF.index])
            # self.allDataDF = pd.concat([self.allDataDF, inputDF/10*6],
                                       # axis=1)

    def getRainfallGPV(self, forecastTime, dateTime=None):
        if dateTime is None:
            dateTime = self.obsDateTime
        grib2_settings = self.getGrib2Params()
        grib2Reader = Grib2Reader(self.inputFilePath)
        grib2FileName = "Z__C_RJTD_" + dateTime.strftime("%Y%m%d%H%M%S") + \
                        "_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"
        # grib2 データから抽出された所望の気象データが格納されたファイル名
        tmpFileName = "tmpfile.csv"

        # grib2 データの情報を表示するためのオプションをリストとして格納します。
        # -V : データの詳細情報を出力するオプションです。
        # -match: そのあとに続く文字列によって切り出し、情報を制約します。
        parameterList = ["-match", "APCP"]

        # grib2 データの情報をオプションを用いて表示します。
        print("******************************************************")
        print(f"Grib2 Data: {grib2FileName}")
        print("a keyward of {}".format(" and ".join(parameterList)))
        print("******************************************************")
        grib2Reader.showInformation(grib2FileName, parameterList)
        gpvDF = pd.DataFrame()
        gpvIndexInterval = 11  # 予測雨量の参照番号の間隔
        initialNum = 21  # 予測雨量の初期参照番号
        endNum = (forecastTime - 1) * gpvIndexInterval + initialNum  # 予測雨量の終了参照番号
        # 参照番号に対応した所望の気象データを１時間ごとに取得するためにループを繰り返します。
        for num in range(initialNum, endNum+1, gpvIndexInterval):
            grib2Num = "1." + str(num)
            # 参照番号を用いて気象データを抽出し、それを csv 形式で保存します。
            print("saving file: {}...".format(tmpFileName))
            grib2Reader.extractData(grib2FileName, grib2Num, tmpFileName)
            # getGPV メソッドの戻り値（リスト）を一時的に変数に格納します。
            tmpDF = grib2Reader.getGPV(tmpFileName,
                                       float(grib2_settings["latitude"]),
                                       float(grib2_settings["longitude"]))
            gpvDF = pd.concat([gpvDF, tmpDF])
            # print(f"data of {tmpList[1]} is detected")
            # dataList.append(tmpList)
        return gpvDF

    def readWaterLevelData(self, used_waterLevelData, timeInterval, timescale,
                           startTime, endTime):
        """
        method for setting data about the water level at the motsukisamu river
        ----
        Input
        ----
        """
        # 望月寒川の時刻水位月報 (1 時間)
        if used_waterLevelData == "motsukisamu" and\
           timeInterval+timescale == "1hours":
            inputFileName = self.inifile.get("inputFile", "waterLevelFileName1")
            sheetList = ["201401", "201402", "201403", "201404", "201405",
                         "201406", "201407", "201408", "201409", "201410",
                         "201411", "201412"]
            inputDF2Dic = pd.read_excel(self.inputFilePath + inputFileName,
                sheet_name=sheetList, skiprows=6, index_col=0)
            startYear = int(startTime.strftime("%Y"))
            startMonth = int(startTime.strftime("%-m"))
            startDay = int(startTime.strftime("%-d"))
            endYear = int(endTime.strftime("%Y"))
            endMonth = int(endTime.strftime("%-m"))
            endDay = int(endTime.strftime("%-d"))
            timestamp = datetime.datetime(startYear, startMonth, startDay)
            readingColumnNo = startDay
            waterLevelList,dateList = [], []
            while timestamp < datetime.datetime(endYear, endMonth, endDay+1):
                # 201409 (2014年9月) の時間水位データをリストに格納
                waterLevelList.append(list(inputDF2Dic[
                    startTime.strftime("%Y%m")].iloc[0:24, readingColumnNo].
                    astype(np.float64)))
                for i in range(24):
                    timestamp += eval("datetime.timedelta(hours=1)")
                    dateList.append(timestamp)
                # timestamp += eval("datetime.timedelta(days=1)")
                readingColumnNo += 1
            waterLevelList = np.array(waterLevelList).flatten()
            tmpSeries = pd.Series(waterLevelList, index=dateList)
            self.allDataDF["water level(obs)"] = tmpSeries[startTime:endTime]

        # 望月寒川の A ロガーデータ (10分)
        elif used_waterLevelData == "motsukisamu" and\
           timeInterval+timescale == "10minutes":
            inputFileName = self.inifile.get("inputFile", "waterLevelFileName2")
            inputDF = pd.read_csv(self.inputFilePath + inputFileName,
                encoding=None, skiprows=None, header=None, sep=",",
                skipinitialspace=True, usecols=[0, 2], index_col=0)
            inputDF.index = pd.to_datetime(inputDF.index)
            inputDF.columns = ["water level(obs)"]
            # self.allDataDF = pd.concat([self.allDataDF, inputDF/1000], axis=1)
            self.allDataDF = pd.concat([self.allDataDF, inputDF/1000], axis=1,
                                       join_axes=[self.allDataDF.index])
