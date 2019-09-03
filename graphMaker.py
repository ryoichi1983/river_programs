# -*- coding: utf-8 -*-
"""
Created on Tue Aug 06 10:35:00 2019

@author: ryoichi.tsurumaki

"""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy as np

class GraphMaker(object):
    """
    class for making and saving the graph
    """
    def __init__(self, allDataDF):
        """
        ----
        Input
        ----
        allDataDF: DataFrame
        """
        self.allDataDF = allDataDF

    def makeFlowRateGraph(self, yNameList, xTickList, xLabel, yLabelList):
        fig, ax = plt.subplots(figsize=(4.86, 3))
        ax1 = ax.twinx()
        ax.grid(True)
        ax1.grid(False)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[0]],
                c="0.1", ls="-", lw=1)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[1]],
                c="0.1", ls="-.", lw=1)
        ax.set_xlim([xTickList[0], xTickList[1]])
        if xTickList[2] == "hours":
            ax.set_ylim([0, 80])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.04, color="0.8")
        elif xTickList[2] == "minutes":
            ax.set_ylim([0, 80])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.01, color="0.8")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.set_xlabel(xLabel, size=9)
        ax.set_ylabel(yLabelList[0], size=9)
        ax1.set_xlabel(xLabel, size=9)
        ax1.set_ylabel(yLabelList[1], size=9)
        ax.tick_params(direction="in", labelsize=9)
        ax1.tick_params(direction="in", labelsize=9)

    def save(self, outputFileName, outputFilePath):
        plt.savefig(outputFilePath + outputFileName, dpi=300)

    def makeWaterLevelGraph(self, yNameList, xTickList, xLabel, yLabelList):
        fig, ax = plt.subplots(figsize=(4.86, 3))
        ax1 = ax.twinx()
        ax.grid(True)
        ax1.grid(False)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[0]],
                c="0.1", ls="-", lw=1)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[1]],
                c="0.1", ls="-.", lw=1)
        ax.set_xlim([xTickList[0], xTickList[1]])
        if xTickList[2] == "hours":
            ax.set_ylim([34, 38])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.04, color="0.8")
        elif xTickList[2] == "minutes":
            ax.set_ylim([34, 38])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.01, color="0.8")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.set_xlabel(xLabel, size=9)
        ax.set_ylabel(yLabelList[0], size=9)
        ax1.set_xlabel(xLabel, size=9)
        ax1.set_ylabel(yLabelList[1], size=9)
        ax.tick_params(direction="in", labelsize=9)
        ax1.tick_params(direction="in", labelsize=9)

    def makeFlow_StorageRelation(self, xName, yName, xLabel, yLabel):
        fig, ax = plt.subplots(figsize=(4.86, 3))
        ax.grid(True)
        ax.plot(self.allDataDF[xName], self.allDataDF[yName],
                c="0.1", marker="o", ms=5, lw=0.5)
        ax.set_xlabel(xLabel, size=9)
        ax.set_ylabel(yLabel, size=9)
        ax.tick_params(direction="in", labelsize=9)

    """
    fig, ax = plt.subplots(figsize=(4.86, 3))
    ax1 = ax.twinx()
    ax.grid(True)
    ax1.grid(False)
    ax.plot(allDataDF.index, allDataDF["water level(obs)"], c="0.1", ls="-", lw=1)
    if cal_settings["timescale"] == "hours":
        ax.set_xlim([datetime.datetime(2014, 9, 10, 0), datetime.datetime(2014, 9, 12, 0)])
        ax.set_ylim([34, 38])
        ax1.set_ylim([100, 0])
        ax1.bar(allDataDF.index, allDataDF["rainfall"], lw=0, width=0.04, color="0.8")
     elif cal_settings["timescale"] == "minutes":
        ax.set_xlim([datetime.datetime(2014, 9, 10, 12), datetime.datetime(2014, 9, 12, 0)])
        ax.set_ylim([34, 38])
        ax1.set_ylim([200, 0])
        ax1.bar(allDataDF.index, allDataDF["rainfall"], lw=0, width=0.01, color="0.8")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlabel("Date", size=9)
    ax.set_ylabel("Water level(m)", size=9)
    ax1.set_xlabel("Date", size=9)
    ax1.set_ylabel("Rainfall (mm/hr)", size=9)
    ax.tick_params(direction="in", labelsize=9)
    ax1.tick_params(direction="in", labelsize=9)
    figureNo = 3
    outputFileName = cal_settings["riverName"] + "_" + \
                     cal_settings["used_flowModel"] + "_" + \
                     cal_settings["used_algorithm"] + "_" + \
                     cal_settings["rainfallStartTime"].strftime("%Y%m%d%H%M") + "_" + \
                     cal_settings["timeInterval"] + \
                     cal_settings["timescale"] + "_" + \
                     str(figureNo) + ".png"
    if cal_settings["used_algorithm"] == "de":
        if cal_settings["used_flowModel"] == "classicOneValueStorageFunc":
            plt.savefig(outputFilePath +
                        "differential_evolution/oneValueStorageFunction/" +
                        outputFileName, dpi=300)
        elif cal_settings["used_flowModel"] == "classicTwoValueStorageFunc":
            plt.savefig(outputFilePath +
                        "differential_evolution/twoValueStorageFunction/" +
                        outputFileName, dpi=300)
        elif cal_settings["used_flowModel"] == "twoStepTwoValueStorageFunc":
            plt.savefig(outputFilePath +
                        "differential_evolution/twoStepTankStorageFunction/" +
                        outputFileName, dpi=300)
        elif cal_settings["used_flowModel"] == "tankModel":
            plt.savefig(outputFilePath +
                        "differential_evolution/tankModel/" +
                        outputFileName, dpi=300)
    elif cal_settings["used_algorithm"] == "kalmanFilter":
        plt.savefig(outputFilePath + "kalmanFilter/" + outputFileName, dpi=300)
    """
