# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:16:00 2019

@author: ryoichi.tsurumaki

"""

import numpy as np
import datetime
from filterpy.kalman import KalmanFilter
from numpy import sqrt
from filterpy.common import Q_discrete_white_noise
import matplotlib.pyplot as plt
from filterpy.common import Q_continuous_white_noise
import statistics
from dataReader import DataReader
from numpy import dot
from scipy.integrate import ode
from storageFuncs import classicTwoValueStorageFunc


class KalmanCalculator:
    """
    カルマンフィルタの計算をおこなうクラスです。
    """
    def __init__(self, cal_settings, catchmentArea):
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
        self.used_objFunc = cal_settings["used_objFunc"]
        self.used_algorithm = cal_settings["used_algorithm"]
        self.used_flowModel = cal_settings["used_flowModel"]
        self.startTime = cal_settings["startTime"]
        self.endTime = cal_settings["endTime"]
        self.timescale = cal_settings["timescale"]
        self.timeInterval = cal_settings["timeInterval"]
        self.forecastTime = cal_settings["forecastTime"]
        self.catchmentArea = catchmentArea
        self.obsName = cal_settings["obsName"]
        self.obsDateTime = cal_settings["obsDateTime"]
        self.prediction = cal_settings["prediction"]
        # self.floodStartTime = cal_settings["floodStartTime"]
        # self.floodEndTime = cal_settings["floodEndTime"]
        # self.rainfallStartTime = cal_settings["rainfallStartTime"]
        # self.rainfallEndTime = cal_settings["rainfallEndTime"]

    def setReadOnlyData(self, readOnlyDataDF, catchmentArea, kalman_settings,
                        baseFlowRate, flowRatio):
        """
        definition of the instance variables
        ----
        Input
        ----
        ----
        """
        self.readOnlyDataDF = readOnlyDataDF
        self.catchmentArea = catchmentArea
        self.baseFlowRate = baseFlowRate
        self.dim_x = kalman_settings["dim_x"]
        self.dim_z = kalman_settings["dim_z"]
        self.dim_u = kalman_settings["dim_u"]
        self.p1 = kalman_settings["p1"]
        self.p2 = kalman_settings["p2"]
        if flowRatio <= 1:
            self.flowRatio = flowRatio
        else:
            self.flowRatio = 1

    def getOutputData(self):
        return self.outputList
    
    # 一段タンク型貯留関数法・拡張カルマンフィルタの場合
    # ----------------------------------------
    def run_EKF(self):
        """
        class filterpy.kalman.KalmanFilter(dim_x, dim_z, dim_u=0)
        Implements a Kalman filter. You are responsible for setting the various 
        state variables to reasonable values; the defaults will not give you a 
        functional filter.
        ----
        Instance Variables
        ----
        dim_x [int] Number of state variables for the Kalman filter. 
              For example, if you are tracking the position and velocity of 
              an object in two dimensions, dim_x would be 4. This is used to
              set the default size of P, Q, and u.
        dim_z [int] Number of of measurement inputs. 
              For example, if the sensor provides you with position in (x,y), 
              dim_z would be 2.
        dim_u [int (optional)] size of the control input, if it is being used. 
              Default value of 0 indicates it is not used.
        compute_log_likelihood [bool (default = True)] Computes log likelihood by default, 
              but this can be a slow computation, so if you never use it you can turn this 
              computation off.
        x [ndarray (dim_x, 1), default = [0,0,0. . . 0]] 
              Current state estimate.
              Any call to update() or predict() updates this variable.
        P [ndarray (dim_x, dim_x), default eye(dim_x)] 
              Current state covariance matrix. 
              Any call to update() or predict() updates this variable.
        x_prior [ndarray (dim_x, 1)] 
              Prior (predicted) state estimate. 
              The *_prior and *_post attributes are for convienence; 
              they store the prior and posterior of the current epoch. Read Only.
        P_prior [ndarray (dim_x, dim_x)] 
              Prior (predicted) state covariance matrix. Read Only.
        x_post [ndarray (dim_x, 1)] 
              Posterior (updated) state estimate. Read Only.
        P_post [ndarray (dim_x, dim_x)] 
              Posterior (updated) state covariance matrix. Read Only.
        z [numpy.array] Last measurement used in update(). Read only.
        R [ndarray (dim_z, dim_z), default eye(dim_x)] measurement uncertainty/noise 
        Q [ndarray (dim_x, dim_x), default eye(dim_x)] Process uncertainty/noise
        F [ndarray (dim_x, dim_x)] state transistion matrix
        H [ndarray (dim_z, dim_x)] measurement function
        B [ndarray (dim_x, dim_u), default 0] optional control transition matrix.
               a value of None will cause the filter to use self.B.
        u [np.array] Optional control vector. If not None, 
               it is multiplied by B to create the control input into the system.
        y [numpy.array] Residual of the update step. Read only.
        K [numpy.array(dim_x, dim_z)] Kalman gain of the update step. Read only.
        S [numpy.array] System uncertainty (P projected to measurement space). Read only.
        log_likelihood [float] log-likelihood of the last measurement.
        likelihood [float] Computed from the log-likelihood. mahalanobis [float] “
        inv [function, default numpy.linalg.inv] 
               If you prefer another inverse function, such as the Moore-Penrose pseudo 
               inverse, set it to that instead: kf.inv = np.linalg.pinv This is only used 
               to invert self.S. If you know it is diagonal, you might choose to set it
               to filterpy.common.inv_diagonal, which is several times faster than numpy.
               linalg.inv for diagonal matrices.
        alpha [float] Fading memory setting.
        z [(dim_z, 1): array_like] measurement for this update. z can be a scalar 
               if dim_z is 1, otherwise it must be convertible to a column vector.
        ----
        Optional Instance Variables
        ----
        alpha : float
        Assign a value > 1.0 to turn this into a fading memory filter.
        ----
        Read-only Instance Variables
        ----
        K [ndarray] Kalman gain that was used in the most recent update() call.
        y [ndarray] Residual calculated in the most recent update() call. 
                    I.e., the different between the measurement and the current 
                    estimated state projected into measurement space (z - Hx)
        S [ndarray] System uncertainty projected into measurement space. 
                    I.e., HPH’ + R. Probably not very useful,
                    but it is here if you want it.
        likelihood [float] Likelihood of last measurment update.
        log_likelihood [float] Log likelihood of last measurment update.
        """
        kalman = KalmanFilter(dim_x=self.dim_x, dim_z=self.dim_z, dim_u=self.dim_u)
        # 降雨開始時刻(初期時刻)
        timestamp = self.startTime
        # 時間グリッド
        gridTime = self.readOnlyDataDF.loc[timestamp, "grid time"]  # timedelta 型
        dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
        # 初期雨量
        rainfallList = [self.readOnlyDataDF.loc[timestamp, "rainfall"] + 10**-10]
        # 初期観測量 (リスト[流量, 雨量])
        # measurableList = [[
        #     self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"]**self.p2 + 10**-10]]
        measurableList = [[
            (self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"] + 10**-10)**self.p2,
            rainfallList[0] + 10**-10]]
        # 状態変数の初期値
        # ----------------------------------------
        kalman.x[0][0] = measurableList[0][0]
        kalman.x[1][0] = 10**-10
        kalman.x[2][0] = 0.9
        kalman.x[3][0] = 1.56
        kalman.x[4][0] = rainfallList[0] + 10**-10
        # kalman.x[4][0] = measurableList[0][1]
        self.kalman_xList = [kalman.x]
        # システム誤差の係数
        sysError = 0.1
        kalman_u = np.array([[1]*self.dim_z]).T
        # kalman_u = 1
        # 状態変数誤差の分散・共分散行列の初期値
        # ----------------------------------------
        kalman.P[0][0] *= (kalman.x[0][0] * sysError)**2
        kalman.P[1][1] *= (kalman.x[1][0] * sysError)**2
        kalman.P[2][2] *= (kalman.x[2][0] * sysError)**2
        kalman.P[3][3] *= (kalman.x[3][0] * sysError)**2
        kalman.P[4][4] *= (kalman.x[4][0] * sysError)**2
        self.kalman_PList = [kalman.P]
        self.kalmanVarianceList = [np.diag(kalman.P)]
        self.outputList = []
        # システム誤差の分散・共分散 (乗算ノイズ)
        # 誤差が小さいとカルマンゲインが小さくなり変数の値がほとんど更新しなくなる。
        kalman_Q = np.array(
            [[0.01, 0, 0, 0, 0],
             [0, 0.01, 0, 0, 0],
             [0, 0, 0.01, 0, 0],
             [0, 0, 0, 0.01, 0],
             [0, 0, 0, 0, 0.01]])

        # カルマンフィルタ計算
        # ----------------------------------------
        while timestamp < self.obsDateTime:
            meanRainfall = statistics.mean(rainfallList)
            lambda1 = 2.8235 * self.catchmentArea**0.24
            lambda2 = 0.2835 * meanRainfall**(-0.2648)
            k11 = lambda1 * self.kalman_xList[-1][3][0]
            k12 = lambda1**2 * lambda2 * self.kalman_xList[-1][3][0]**2
            kalman_F = kalman.FJacobian(dt, k11, k12, self.p1, self.p2)
            kalman_B = kalman.calcControlVector(dt, k11, k12, self.p1, self.p2)
            # 観測の相対誤差 (10% 仮定：最小二乗法によって誤差を求める必要あり)
            # obsError = sqrt(self.p2 * (self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"] +
            #                            10**-10)**(self.p2 - 1) * 1.1)
            obsError = 0.1
            # 観測行列
            # kalman_H = np.array([[1, 0, 0, 0, 0]])
            kalman_H = np.array([[1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])
            kalman.predict(B=kalman_B, u=kalman_u, F=kalman_F, Q=kalman_Q)
            timestamp += gridTime
            gridTime = self.readOnlyDataDF.loc[timestamp, "grid time"]  # timedelta 型
            dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
            rainfallList.append(self.readOnlyDataDF.loc[timestamp, "rainfall"])
            # measurableList.append([
            #         self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"]**self.p2 +
            #         10**-10])
            measurableList.append([
                    self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"]**self.p2 +
                    10**-10, rainfallList[-1] + 10**-10])
            # 観測ノイズの分散・共分散
            # kalman_R = np.array([[obsError * self.kalman_xList[-1][0][0]]])
            kalman_R = np.array([[obsError * self.kalman_xList[-1][0][0], 0],
                                 [0, obsError * self.kalman_xList[-1][4][0]]])
            kalman.update(z=measurableList[-1], R=kalman_R, H=kalman_H)  # 更新(z=[[None]])
            self.kalman_xList.append(kalman.x)
            self.kalman_PList.append(kalman.P)
            self.kalmanVarianceList.append(np.diag(kalman.P))
            # システム誤差の分散・共分散 (乗算ノイズ)
            # 雨量の場合は予測誤差とする
            kalman_Q = np.array(
                    [[(self.kalman_xList[-1][0][0] * sysError)**2, 0, 0, 0, 0],
                     [0, (self.kalman_xList[-1][1][0] * sysError)**2, 0, 0, 0],
                     [0, 0, (self.kalman_xList[-1][2][0] * sysError)**2, 0, 0],
                     [0, 0, 0, (self.kalman_xList[-1][3][0] * sysError)**2, 0],
                     [0, 0, 0, 0, (self.kalman_xList[-1][4][0] * sysError)**2]])
            flowRate = self.kalman_xList[-1][0][0]**(1/self.p2)
            errorFlowRate = self.catchmentArea / 3.6 / self.p2 * \
                self.kalman_xList[-1][0][0]**(1/self.p2 - 1) * \
                self.kalman_PList[0][0][0]**0.5 + self.baseFlowRate
            # sigma_WL = sigma_flow / (4 * flowRate)  要修正
            storageHeight = k11 * flowRate**self.p1 + k12 * self.kalman_xList[-1][1][0]
            self.outputList.append([timestamp, flowRate, rainfallList[-1], 
                                    errorFlowRate, storageHeight])
        # time = [i for i in range(len(kalman_xList))]
        # kalman_xList = np.array(kalman_xList).reshape(25, 5)
        # plt.plot(time, kalman_xList)
        if self.prediction:
            dataReader_kal = DataReader()
            dataReader_kal.setInputFilePath()
            grib2_settings = dataReader_kal.getGrib2Params()
            forecastDateTime = self.obsDateTime + \
                eval("datetime.timedelta(" + grib2_settings["timescale"] +
                     "=" + self.forecastTime + ")")
            gridTime = eval("datetime.timedelta(" + grib2_settings["timescale"] + "=" +
                        grib2_settings["timeInterval"] + ")")
            rainfallGPVDF = dataReader_kal.getRainfallGPV(
                forecastTime=int(self.forecastTime), dateTime=self.obsDateTime)
            # rainfallList = []
            # for rainfall in rainfallGPVDF["rainfall"].values.tolist():
                # rainfallList.append(rainfall)
            self.forecast(kalman, forecastDateTime, gridTime, rainfallGPVDF)
        """
        if int(self.forecastTime) != 0:
            dataReader_kal = DataReader()
            dataReader_kal.setInputFilePath()
            grib2_settings = dataReader_kal.getGrib2Params()
            forecastDateTime = self.obsDateTime + \
                eval("datetime.timedelta(" + grib2_settings["timescale"] +
                     "=" + self.forecastTime + ")")
            gridTime = eval("datetime.timedelta(" + grib2_settings["timescale"] + "=" +
                        grib2_settings["timeInterval"] + ")")
            rainfallGPVDF = dataReader_kal.getRainfallGPV(
                forecastTime=int(self.forecastTime), dateTime=self.obsDateTime)
            # rainfallList = []
            # for rainfall in rainfallGPVDF["rainfall"].values.tolist():
                # rainfallList.append(rainfall)
            self.forecast(kalman, forecastDateTime, gridTime, rainfallGPVDF)
        """

    def forecast(self, kalman, forecastDateTime, gridTime, rainfallGPVDF=None):
        """
        class for the forecast (15hours) of state vector x and error function variance P
        ----
        Input
        ----
        """
        if rainfallGPVDF is None:
            rainfallGPVDF = self.rainfallGPVDF
        meanRainfall = (self.readOnlyDataDF.loc[self.startTime:self.obsDateTime,
                                                "rainfall"].sum() + \
               rainfallGPVDF["rainfall"].sum()) / \
               (len(self.readOnlyDataDF.loc[self.startTime:self.obsDateTime,
                                            "rainfall"]) + \
                len(rainfallGPVDF["rainfall"]))
        dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
        sol = ode(classicTwoValueStorageFunc)
        sol.set_integrator("dop853")
        t0 = 0
        y0 = [self.kalman_xList[-1][0][0], self.kalman_xList[-1][1][0]]
        flowRatio = self.kalman_xList[-1][2][0]
        lambda1 = 2.8235 * self.catchmentArea**0.24
        lambda2 = 0.2835 * meanRainfall**(-0.2648)
        k11 = lambda1 * self.kalman_xList[-1][3][0]
        k12 = lambda1**2 * lambda2 * self.kalman_xList[-1][3][0]**2
        sol.set_initial_value(y=y0, t=t0)
        timestamp = self.obsDateTime + gridTime
        # sysError = 0.1
        # kalman_u = np.array([[1]*self.dim_z]).T
        while timestamp <= forecastDateTime:
            rainfall = rainfallGPVDF.loc[timestamp, "rainfall"] * flowRatio
            sol.set_f_params(k11, k12, self.p1, self.p2, rainfall)
            integrateResultList = sol.integrate(sol.t+dt)
            flowRate = (integrateResultList[0])**(1/self.p2) * self.catchmentArea / 3.6 + \
                       self.baseFlowRate
            storageHeight = k11 * integrateResultList[0]**(self.p1/self.p2) + \
                            k12 * integrateResultList[1]
            errorFlowRate = np.nan
            self.outputList.append([timestamp, flowRate, 
                                    rainfallGPVDF.loc[timestamp, "rainfall"],
                                    errorFlowRate, storageHeight])
            timestamp += gridTime

    """
    def forecast(self, kalman, rainfallList=None, obsDateTime=None, forecastTime=None):
        class for the forecast (15hours) of state vector x and error function variance P
        # ----
        # Input
        # ----
        if obsDateTime is None:
            obsDateTime = self.obsDateTime
        if forecastTime is None:
            forecastTime = self.forecastTime
        if rainfallList is None:
            rainfallList = []
        dataReader_kal = DataReader()
        dataReader_kal.setInputFilePath()
        grib2_settings = dataReader_kal.getGrib2Params()
        rainfallGPVDF = dataReader_kal.getRainfallGPV(
            forecastTime=int(self.forecastTime), dateTime=obsDateTime)
        for rainfall in rainfallGPVDF["rainfall"].values.tolist():
            rainfallList.append(rainfall)
        meanRainfall = statistics.mean(rainfallList)
        forecastDateTime = obsDateTime + \
            eval("datetime.timedelta(" + grib2_settings["timescale"] +
                 "=" + self.forecastTime + ")")
        gridTime = eval("datetime.timedelta(" + grib2_settings["timescale"] + "=" +
                        grib2_settings["timeInterval"] + ")")
        dt = gridTime.days * 24 + gridTime.seconds / 3600  # int 型
        kalman_u = np.array([[1]*self.dim_z]).T
        sysError = 0.1
        timestamp = obsDateTime + gridTime
        while timestamp <= forecastDateTime:
            # meanRainfall = statistics.mean(rainfallList)
            # meanRainfall = rainfallGPVDF["rainfall"].mean() + 10**-10
            lambda1 = 2.8235 * self.catchmentArea**0.24
            lambda2 = 0.2835 * meanRainfall**(-0.2648)
            k11 = lambda1 * self.kalman_xList[-1][3][0]
            k12 = lambda1**2 * lambda2 * self.kalman_xList[-1][3][0]**2
            kalman.x[4][0] = rainfallGPVDF.loc[timestamp, "rainfall"]
            kalman_F = kalman.FJacobian(dt, k11, k12, self.p1, self.p2)
            kalman_B = kalman.calcControlVector(dt, k11, k12, self.p1, self.p2)
            # 観測の相対誤差 (10% 仮定：最小二乗法によって誤差を求める必要あり)
            # obsError = sqrt(self.p2 * self.readOnlyDataDF.loc[timestamp, "flow rate(HQ)"]**(
            #     self.p2 - 1) * 1.1)
            obsError = 0.1
            # システム誤差の分散・共分散 (乗算ノイズ)
            # 雨量の場合は予測誤差とする
            kalman_Q = np.array(
                    [[(self.kalman_xList[-1][0][0] * sysError)**2, 0, 0, 0, 0],
                     [0, (self.kalman_xList[-1][1][0] * sysError)**2, 0, 0, 0],
                     [0, 0, (self.kalman_xList[-1][2][0] * sysError)**2, 0, 0],
                     [0, 0, 0, (self.kalman_xList[-1][3][0] * sysError)**2, 0],
                     [0, 0, 0, 0, (self.kalman_xList[-1][4][0] * sysError)**2]])
            # print((dot(kalman_F, kalman.x) + dot(kalman_B, kalman_u))[0][0])
            kalman.predict(B=kalman_B, u=kalman_u, F=kalman_F, Q=kalman_Q)
            self.kalman_xList.append(kalman.x)
            self.kalman_PList.append(kalman.P)
            self.kalmanVarianceList.append(np.diag(kalman.P))
            flowRate = self.kalman_xList[-1][0][0]**(1/self.p2)
            errorFlowRate = self.catchmentArea / 3.6 / self.p2 * \
                self.kalman_xList[-1][0][0]**(1/self.p2 - 1) * \
                self.kalman_PList[0][0][0]**0.5
            # sigma_WL = sigma_flow / (4 * flowRate)  要修正
            self.outputList.append([timestamp, flowRate, 
                                    rainfallGPVDF.loc[timestamp, "rainfall"],
                                    errorFlowRate])
            timestamp += gridTime
    """
