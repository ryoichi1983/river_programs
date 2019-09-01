# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:41:00 2019

@author: ryoichi.tsurumaki
"""

import subprocess
import os
import pandas as pd


class Grib2Reader:
    """
    grib2 (bin) ファイルを wgrib2 モジュールでデコードし、所望の気象データを
    読み込むクラスを定義します。
    """

    def __init__(self, inputDirPath, outputDirPath="./tmpDir/"):
        """
        ----
        Input
        ----
        * inputDirPath: str
            input file directory
        * outputDirPath: str
            save directory path for output file
        ----
        """
        self.inputDirPath = inputDirPath
        self.outputDirPath = outputDirPath
        print("#########################")
        print("output directory")
        print(outputDirPath)
        print("#########################")
        # 結果のダウンロード先のディレクトリを作成します。
        subprocess.run(['mkdir', '-p', self.outputDirPath])

    def showInformation(self, grib2FileName, parameterList):
        """
        grib2 ファイルに存在する気象データの情報を表示するメソッドを定義します。
        ----
        Input
        ----
        * grib2FileName: str
            grib2 file name
        * parameterList: list
            wgrib2 command options
        ----
        """

        # grib2 気象データの情報を表示するためのコマンドをリストで定義します。
        commandList = parameterList[:]
        commandList.insert(0, "wgrib2")
        commandList.append(self.inputDirPath + grib2FileName)

        # grib2 気象データの情報を表示します。
        subprocess.run(commandList)

    def extractData(self, grib2FileName, grib2Num, tmpFileName, 
                    outputFileType="csv"):
        """
        grib2 ファイルから所望の気象データが保存されているファイルを csv 形式で
        抽出します。
        ----
        Input
        ----
        * grib2FileName: str
            grib2 file name for input
        * grib2Num: str
            reference number of data for grib2 file
        * outputFileType: str
            csv or txt
        * outputFileName: str
            csv or txt
        ----
        """

        # grib2 気象データの情報を表示するためのコマンドをリストで定義します。
        commandList = ["wgrib2", self.inputDirPath + grib2FileName, "-d",
                       grib2Num, "-" + outputFileType, self.outputDirPath +
                       tmpFileName]
        # grib2 気象データを保存します。
        subprocess.run(commandList)

    def getGPV(self, tmpFileName, latitude, longitude):
        """
        grib2 ファイルから抽出した csv ファイルから所望のデータを取り出し、
        日時や物理量などが数値文字列として格納されたリストを返します。
        ----
        Input
        ----
        * tmpFileName: str
            extracted meteolorogical data file name by function extractData
        ----
        """

        # grib2 データから抽出した csv ファイルの File オブジェクトを読み込み
        # モードで作成します。
        """
        with open(self.outputDirPath + tmpFileName, "r", encoding="utf8") \
        as fileObj:
            # 行単位で文字列要素を取り出し、１行ずつ文末までループを繰り返します。
            for i, line in enumerate(fileObj):
                # 文字列要素が改行のみの行はスキップします。
                if line == "\n":
                    continue
                # 文字列要素の行末の改行コードを取り除きます。
                line = line.rstrip()
                # カンマ区切りの文字列要素をリストに変換します。
                lineList = line.split(",")
                # 所望の緯度・経度に等しい場合
                if (float(lineList[4]) == longitude) and \
                   (float(lineList[5]) == latitude):
                    return lineList
        """
        outputDF = pd.read_csv(self.outputDirPath + tmpFileName,
                               encoding="utf-8", skiprows=None, header=None,
                               sep=",", skipinitialspace=True,
                               usecols=[1, 4, 5, 6], index_col=0)
        outputDF.columns = ["longitude", "latitude", "rainfall"]
        outputDF.index = pd.to_datetime(outputDF.index)
        # gpv = outputDF.query(eval("longitude==" + longitude + " & latitude==" + latitude))
        # gpv = outputDF.query("longitude==141.375 & latitude==43.05")
        gpv = outputDF[(outputDF['longitude'] == longitude) &
                       (outputDF['latitude'] == latitude)]
        return gpv

# ----------------------------------------------------------------------------
# main 関数の定義
# ----------------------------------------------------------------------------
if __name__ == '__main__':

    # grib2 ファイルが保存されているディレクトリ名
    inputDirPath = \
    "./inputFiles/"

    # grib2 ファイル名
    grib2FileName = \
    "Z__C_RJTD_20140910000000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"

    # 結果を保存するディレクトリ名
    outputDirPath = \
    "./outputFiles/"

    # Grib2Reader 型のインスタンスの作成
    grib2Reader = Grib2Reader(inputDirPath, outputDirPath)

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

    latitude = 43.05  # 流域緯度
    longitude = 141.375  # 流域経度

    # 気象データを格納するためのリスト変数を定義します。
    dataList = []
    
    # 参照番号に対応した所望の気象データを1時間ごとに15時間分取得するた
    # めにループを繰り返します。
    gpvDF = pd.DataFrame()
    for num in range(21, 176, 11):
        grib2Num = "1." + str(num)

        # 参照番号を用いて気象データを抽出し、それを csv 形式で保存します。
        print("saving file: {}...".format(tmpFileName))
        grib2Reader.extractData(grib2FileName, grib2Num, tmpFileName)

        # getGPV メソッドの戻り値（リスト）を一時的に変数に格納します。
        tmpDF = grib2Reader.getGPV(tmpFileName, latitude, longitude)
        gpvDF = pd.concat([gpvDF, tmpDF])

        # print(f"data of {tmpList[1]} is detected")
        # dataList.append(tmpList)
    print(dataList)
    
