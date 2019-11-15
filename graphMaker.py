# -*- coding: utf-8 -*-
"""
Created on Tue Aug 06 10:35:00 2019

@author: ryoichi.tsurumaki

"""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import random
import math
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

class GraphMaker(object):
    """
    class for making and saving the graph
    """
    def __init__(self, allDataDF=None):
        """
        ----
        Input
        ----
        allDataDF: DataFrame
        """
        if allDataDF is not None:
            self.allDataDF = allDataDF
        else:
            self.allDataDF = []

    def makeFlowRateGraph(self, yNameList, xTickList, xLabel, yLabelList):
        fig, ax = plt.subplots(figsize=(4.86, 3))
        ax1 = ax.twinx()
        ax.grid(True)
        ax1.grid(False)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[0]],
                c="0.1", ls="-", lw=1)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[1]],
                c="0.1", ls="--", lw=1)
        ax.set_xlim([xTickList[0], xTickList[1]])
        if xTickList[2] == "hours":
            ax.set_ylim([0, 70])
            ax1.set_ylim([140, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.04, color="0.8")
            ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
        elif xTickList[2] == "minutes":
            ax.set_ylim([0, 80])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=0.01, color="0.8")
            ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
        elif xTickList[2] == "days":
            ax.set_ylim([0, 30])
            ax1.set_ylim([200, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[2]],
                    lw=0, width=1, color="0.8")
            ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
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
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=24))
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

    def makeAnimation(self, populationList, population_energyList):
        fig =  plt.figure(figsize=(4.86, 3))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel("x", size = 9)
        ax.set_ylabel("y", size = 9)
        ax.set_zlabel("z", size = 9)
        ax.set_zlim([7, 8.8])
        # ax.set_xlim([-1.5, 1.5])
        # ax.set_ylim([-1.5, 1.5])

        images = []
        for populations, z in zip(populationList, population_energyList):
            x, y = [], []
            for xy in populations:
                x.append(xy[0])
                y.append(xy[1])
            # for xy, z in zip(populations, population_energies):
                # image = ax.scatter(xy[0], xy[1], z)
                # images.append([image])
            x = np.array(x)
            y = np.array(y)
            z = z.flatten()
            image = ax.scatter(x, y, z)
            images.append([image])
        ani = animation.ArtistAnimation(fig, 
                                        images, 
                                        blit=True, 
                                        repeat=False,
                                        interval=500)
        ani.save("./results/output.gif", writer="imagemagick")
        plt.show()

    def makeWLRFGraph(self, yNameList, xTickList, xLabel, yLabelList):
        fig, ax = plt.subplots(figsize=(4.86, 3))
        ax1 = ax.twinx()
        ax.grid(True)
        ax1.grid(False)
        ax.plot(self.allDataDF.index, self.allDataDF[yNameList[0]],
                c="0.1", ls="-", lw=1)
        ax.set_xlim([xTickList[0], xTickList[1]])
        if xTickList[2] == "hours":
            ax.set_ylim([34, 38])
            ax1.set_ylim([80, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[1]],
                    lw=0, width=0.04, color="0.8")
        elif xTickList[2] == "minutes":
            ax.set_ylim([34, 38])
            ax1.set_ylim([160, 0])
            ax1.bar(self.allDataDF.index, self.allDataDF[yNameList[1]],
                    lw=0, width=0.01, color="0.8")
        # ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=6))
        ax.set_xlabel(xLabel, size=9)
        ax.set_ylabel(yLabelList[0], size=9)
        ax1.set_xlabel(xLabel, size=9)
        ax1.set_ylabel(yLabelList[1], size=9)
        ax.tick_params(direction="in", labelsize=9)
        ax1.tick_params(direction="in", labelsize=9)


if __name__ == '__main__':
    """
    # 更新する内容
    def _update_plot(i, fig, im):
        rad = math.radians(i)

        # 前回のフレーム内容を一旦削除
        if len(im) > 0:
            im[0].remove()
            im.pop()

        im.append(plt.scatter(math.cos(rad), math.sin(rad)))
    fig =  plt.figure()

    # グラフを中央に表示
    ax = fig.add_subplot(1,1,1)

    # グラフの目盛範囲設定
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])

    im = [] # フレーム更新の際に前回のプロットを削除するために用意

    """
    graphMaker = GraphMaker()
    graphMaker.makeAnimation()

    # 表示
    #plt.show()

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
