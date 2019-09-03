# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:16:00 2019

@author: ryoichi.tsurumaki

"""

import numpy as np
from objectivefunctions import rmse
from objectivefunctions import nashsutcliffe
from tankModelCalculator import TankModelCalculator
from storageFuncs import classicOneValueStorageFunc
from storageFuncs import classicTwoValueStorageFunc
from storageFuncs import twoStepTwoValueStorageFunc
from scipy.integrate import ode
from filterpy.common import Q_discrete_white_noise
from filterpy.common import Q_continuous_white_noise


class ParamOptimizer:
    """
    パラメータ同定のための最適化をおこなうクラスです。
    """
    def __init__(self, calc_settings):
        """
        ----
        Input
        ----
        used_objFunc: str
            * objective function for searching for parameters
            * rmse (root mean square)
            * nash (nash-sutcliffe efficiency)
            * chi (chi squared test)
        used_model: str
            * flow model for the river
            * tankModel or classicOneValueStorageFunc or classicTwoValueStorageFunc
              or twoStepTwoValueStorageFunc
        used_algorithm: str
            * optimization model to find a minimization of function
            * de(differential evolution) or sceua(SCE-UA method)
        """
        self.used_objFunc = calc_settings["used_objFunc"]
        self.used_algorithm = calc_settings["used_algorithm"]
        self.used_flowModel = calc_settings["used_flowModel"]
        self.startTime = calc_settings["startTime"]
        self.endTime = calc_settings["endTime"]
        self.timescale = calc_settings["timescale"]
        self.timeInterval = calc_settings["timeInterval"]

    def calculateObjfunction(self, simulation, evaluation):
        """
        objective fuction for the optimization
        ----
        Input
        ----
        simulation: list
            simulation data to compared with evaluation data
        evaluation: list
            Observed data to compared with simulation data.
        objectivefunction: float
            calculated value by an objective function
        ----
        """
        if self.used_objFunc == "rmse":
            objectivefunction = rmse(evaluation=evaluation, simulation=simulation)
        elif self.used_objFunc == "nash":
            objectivefunction = nashsutcliffe(evaluation=evaluation, simulation=simulation)
        elif self.used_objFunc == "chi":
            pass
        return objectivefunction

    def setReadDataForTank(self, readDataDF, catchmentArea, tankModel_settings):
        """
        definition of the instance variables
        ----
        Input
        ----
        readDataDF: DataFrame
            input data for the flow analysis (READ ONLY!!)
        catchmentArea: float
            area to catch the rain for the river (km^2)
        tankModel_settings: dic(float)
            parameters for tank model
        ----
        """
        self.readDataDF = readDataDF
        self.catchmentArea = catchmentArea
        self.tankModel_settings = tankModel_settings

    def setTankParams(self, x):
        """
        ----
        Input
        ----
        x: list
        alpha1: 第 1 タンクの下部横孔の流出係数(hr^-1)
        alpha2: 第 1 タンクの上部横孔の流出係数(hr^-1)
        alpha3: 第 2 タンクの横孔の流出係数(hr^-1)
        alpha4: 第 3 タンクの横孔の流出係数(hr^-1)
        beta1: 第 1 タンクの下孔の浸透係数(hr^-1)
        beta2: 第 1 タンクの下孔の浸透係数(hr^-1)
        beta3: 第 2 タンクの下孔の浸透係数(hr^-1)
        height1: 第 1 タンクの下部横孔の高さ(mm)
        height2: 第 1 タンクの上部横孔の高さ(mm)
        height3: 第 2 タンクの横孔の高さ(mm)
        height4: 第 3 タンクの横孔の高さ(mm)
        storage1: 第 1 タンクの初期貯留高(mm)
        storage2: 第 2 タンクの初期貯留高(mm)
        storage3: 第 3 タンクの初期貯留高(mm)
        ----
        """
        self.tankModel_settings["alpha1"] = x[0]
        self.tankModel_settings["alpha2"] = x[1]
        self.tankModel_settings["alpha3"] = x[2]
        self.tankModel_settings["alpha4"] = x[3]
        self.tankModel_settings["beta1"] = x[4]
        self.tankModel_settings["beta2"] = x[5]
        self.tankModel_settings["beta3"] = x[6]
        self.tankModel_settings["height1"] = x[7]
        self.tankModel_settings["height2"] = x[8]
        self.tankModel_settings["height3"] = x[9]
        self.tankModel_settings["height4"] = x[10]
        self.tankModel_settings["storage1"] = 0
        self.tankModel_settings["storage2"] = 0
        self.tankModel_settings["storage3"] = 0

    """
    def getTankParams(self, x):
        self.setTankParams(x)
        return self.tankModel_settings
    """

    def getOutputData(self):
        return self.outputList
    
    def searchFunc(self, x):
        """
        function for the global minimization
        ----
        Input
        ----
        x: list
            parameters for optimizing the function
        ----
        return: numerical value of the objective function
        ----
        """
        simulationList = []
        self.outputList = []
        # タンクモデルの場合
        # ----------------------------------------
        if self.used_flowModel == "tankModel":
            evaluationList = []
            self.setTankParams(x)
            tankModelCalculator = TankModelCalculator(**self.tankModel_settings)
            timestamp = self.startTime
            while timestamp <= self.endTime:
                gridTime = self.readDataDF.loc[timestamp, "grid time"]  # timedelta 型
                dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
                rainfall = self.readDataDF.loc[timestamp, "rainfall"]
                self.outputList.append(tankModelCalculator.calculateTank(dt, rainfall))
                flowRate = self.outputList[-1][4] * self.catchmentArea / 3.6
                self.outputList[-1] = np.append(self.outputList[-1],
                                                [timestamp, flowRate, rainfall])
                simulationList.append(flowRate)
                evaluationList.append(self.readDataDF.loc[timestamp, "flow rate(HQ)"])
                timestamp += gridTime
            return self.calculateObjfunction(simulationList, evaluationList)

        # 古典的一価非線形貯留関数の場合
        # ----------------------------------------
        elif self.used_flowModel == "classicOneValueStorageFunc":
            evaluationList = []
            sol = ode(classicOneValueStorageFunc)
            sol.set_integrator("dop853")
            y0, t0 = 0, 0
            sol.set_initial_value(y=y0, t=t0)
            # timestamp = self.rainfallStartTime
            timestamp = self.startTime
            a, p = x[0], x[1]
            # while sol.successful() and timestamp <= self.endTime:
            while timestamp <= self.endTime:
                gridTime = self.readDataDF.loc[timestamp, "grid time"]
                dt = gridTime.days * 24 + gridTime.seconds / 3600
                rainfall = self.readDataDF.loc[timestamp, "effective rainfall"]
                sol.set_f_params(a, p, rainfall)
                integrateResultList = sol.integrate(sol.t+dt)
                flowRate = (integrateResultList[0])**(1/p) * self.catchmentArea / 3.6
                storageHeight = a * integrateResultList[0]
                if storageHeight < 0:
                    return np.inf
                self.outputList.append([timestamp, flowRate, rainfall, storageHeight])
                simulationList.append(flowRate)
                evaluationList.append(self.readDataDF.loc[timestamp, "flow rate(HQ)"])
                timestamp += gridTime
            return self.calculateObjfunction(simulationList, evaluationList)
            """
            if sol.successful():
                return self.calculateObjfunction(simulationList, evaluationList)
            else:
                print("## INTEGRATION ERROR ##")
                # return "error"
                return 10**20  # 誤差の上限値(根拠なし)
            """

        # 古典的二価非線形貯留関数の場合
        # ----------------------------------------
        elif self.used_flowModel == "classicTwoValueStorageFunc":
            evaluationList = []
            sol = ode(classicTwoValueStorageFunc)
            sol.set_integrator("dop853")
            y0 = [0, 0]  # (x, y)
            t0 = 0
            sol.set_initial_value(y=y0, t=t0)
            timestamp = self.startTime
            a, b, m, n = x[0], x[1], x[2], x[3]
            while timestamp <= self.endTime:
            # while sol.successful() and timestamp <= self.endTime:
                gridTime = self.readDataDF.loc[timestamp, "grid time"]  # timedelta 型
                dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型(時間)
                rainfall = self.readDataDF.loc[timestamp, "effective rainfall"]
                sol.set_f_params(a, b, m, n, rainfall)
                integrateResultList = sol.integrate(sol.t+dt)
                flowRate = (integrateResultList[0])**(1/n) * self.catchmentArea / 3.6
                storageHeight = a * integrateResultList[0]**(m/n) + \
                                b * integrateResultList[1]
                if storageHeight < 0:
                    return np.inf
                # if not sol.successful():
                    # return np.inf
                self.outputList.append([timestamp, flowRate, rainfall, storageHeight])
                simulationList.append(flowRate)
                evaluationList.append(self.readDataDF.loc[timestamp, "flow rate(HQ)"])
                timestamp += gridTime
            return self.calculateObjfunction(simulationList, evaluationList)
            """
            if sol.successful():
                return self.calculateObjfunction(simulationList, evaluationList)
            else:
                print("## INTEGRATION ERROR ##")
                return np.inf
            """

        # 二段タンク型二価非線形貯留関数の場合
        # ----------------------------------------
        elif self.used_flowModel == "twoStepTwoValueStorageFunc":
            evaluationList = []
            # evaluationList = list(self.readDataDF["flow rate(HQ)"].loc[
                # self.rainfallStartTime: self.floodEndTime])
            sol = ode(twoStepTwoValueStorageFunc)
            sol.set_integrator("dop853")
            y0 = [0, 0, 0, 0]  # (y1, y2, y3, y4)
            t0 = 0
            sol.set_initial_value(y=y0, t=t0)
            timestamp = self.startTime
            k11 = x[0]
            k12 = x[1]
            k13 = x[2]
            k21 = x[3]
            k22 = x[4]
            p1 = x[5]
            p2 = x[6]
            while timestamp <= self.endTime:
            # while sol.successful() and timestamp <= self.endTime:
                gridTime = self.readDataDF.loc[timestamp, "grid time"]  # timedelta 型
                dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
                # gridTime = self.readDataDF.loc[timestamp, "grid time"]
                rainfall = self.readDataDF.loc[timestamp, "rainfall"]
                sol.set_f_params(k11, k12, k13, k21, k22, p1, p2, rainfall)
                integrateResultList = sol.integrate(sol.t+dt)
                surfaceToGround = integrateResultList[0]**(1/p2) /\
                                  integrateResultList[2]
                # 差分計算の出力 (リストNo.0 が求めたい量であることに注意)
                flowRate = (integrateResultList[0]**(1/p2) + integrateResultList[2]) * \
                    self.catchmentArea / 3.6
                storageHeight = k11 * integrateResultList[0]**(p1/p2) + \
                                k12 * integrateResultList[1] + \
                                k21 * integrateResultList[2] + \
                                k22 * integrateResultList[3]
                if storageHeight < 0:
                    return np.inf
                self.outputList.append([timestamp, flowRate, rainfall, storageHeight])
                simulationList.append(flowRate)
                evaluationList.append(self.readDataDF.loc[timestamp, "flow rate(HQ)"])
                timestamp += gridTime
            return self.calculateObjfunction(simulationList, evaluationList)
            """
            if sol.successful():
                return self.calculateObjfunction(simulationList, evaluationList)
            else:
                print("## INTEGRATION ERROR ##")
                return 10**10  # 誤差の上限値(根拠なし)
            """

    def getBounds(self):
        """
        definition of parameter ranges for the optimization
        ----
        Input
        ----
        ----
        """
        if self.used_flowModel == "tankModel":
            return [(0, 10), (0, 10), (0, 10), (0, 10), (0, 10), (0, 10), 
                    (0, 10), (0, 1000), (0, 1000), (0, 1000), (0, 1000)]
        elif self.used_flowModel == "classicOneValueStorageFunc":
            return [(0, 10), (0, 1.5)]
        elif self.used_flowModel == "classicTwoValueStorageFunc":
            return [(0, 10), (0, 10), (0, 1), (0, 1)]
        elif self.used_flowModel == "twoStepTwoValueStorageFunc":
            return [(0, 10), (0, 10), (0, 10), (0, 10), (0, 10), (0, 1), (0, 1)]

