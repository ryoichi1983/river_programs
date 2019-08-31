# -*- coding: utf-8 -*-
"""
Created on Tue Aug 06 10:35:00 2019

@author: ryoichi.tsurumaki

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
import configparser
from tankModelCalculator import TankModelCalculator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from paramOptimizer import ParamOptimizer
from scipy.optimize import differential_evolution
from dataReader import DataReader


class GraphMaker(object):
    """
    class for making graph
    """
    def __init__(self, allDataDF):
        self.allDataDF = allDataDF

    fig, ax = plt.subplots(figsize=(4.86, 3))
    ax1 = ax.twinx()
    ax.grid(True)
    ax1.grid(False)
    ax.plot(allDataDF.index, allDataDF["flow rate(HQ)"], c="0.1", ls="-", lw=1)
    ax.plot(allDataDF.index, allDataDF["flow rate(cal)"], c="0.1", ls="-.", lw=1)
    if cal_settings["timescale"] == "hours":
#        ax.set_xlim([datetime.datetime(2014, 9, 10, 0), datetime.datetime(2014, 9, 12, 0)])
        ax.set_xlim([cal_settings["startTime"], cal_settings["endTime"]])
        ax.set_ylim([0, 50])
        ax1.set_ylim([100, 0])
        ax1.bar(allDataDF.index, allDataDF["rainfall"], lw=0, width=0.04, color="0.8")
    elif cal_settings["timescale"] == "minutes":
        ax.set_xlim(cal_settings["startTime"], cal_settings["endTime"])
#        ax.set_xlim([datetime.datetime(2014, 9, 10, 12), datetime.datetime(2014, 9, 12, 0)])
        ax.set_ylim([0, 80])
        ax1.set_ylim([200, 0])
        ax1.bar(allDataDF.index, allDataDF["rainfall"], lw=0, width=0.01, color="0.8")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d-%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlabel("Date", size=9)
    ax.set_ylabel("flow rate (m$^3$/s)", size=9)
    ax1.set_xlabel("Date", size=9)
    ax1.set_ylabel("Rainfall (mm/hr)", size=9)
    ax.tick_params(direction="in", labelsize=9)
    ax1.tick_params(direction="in", labelsize=9)
    figureNo = 1
    outputFileName = cal_settings["riverName"] + "_" + \
                     cal_settings["used_flowModel"] + "_" + \
                     cal_settings["used_algorithm"] + "_" + \
                     cal_settings["rainfallStartTime"].strftime("%Y%m%d%H%M") + "_" + \
                     cal_settings["timeInterval"] + \
                     cal_settings["timescale"] + "_" + \
                     str(figureNo) + ".png"
    """
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

    fig, ax = plt.subplots(figsize=(4.86, 3))
    ax1 = ax.twinx()
    ax.grid(True)
    ax1.grid(False)
    ax.plot(allDataDF.index, allDataDF["water level(obs)"], c="0.1", ls="-", lw=1)
    ax.plot(allDataDF.index, allDataDF["water level(cal)"], c="0.1", ls="-.", lw=1)
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
    figureNo = 2
    outputFileName = cal_settings["riverName"] + "_" + \
                     cal_settings["used_flowModel"] + "_" + \
                     cal_settings["used_algorithm"] + "_" + \
                     cal_settings["rainfallStartTime"].strftime("%Y%m%d%H%M") + "_" + \
                     cal_settings["timeInterval"] + \
                     cal_settings["timescale"] + "_" + \
                     str(figureNo) + ".png"
    """
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

    """
    fig, ax = plt.subplots(figsize=(4.86, 3))
    ax.grid(True)
    # ax.scatter(allDataDF["flow rate(cal)"], allDataDF["storage height"], c="0.1", s=10)
    ax.plot(allDataDF["flow rate(cal)"], allDataDF["storage height"], c="0.1", marker="o", ms=5, lw=0.5)
    ax.set_xlabel("flow rate(m$^3$/s)", size=9)
    ax.set_ylabel("storage height(mm)", size=9)
    ax.tick_params(direction="in", labelsize=9)
    figureNo = 4
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

