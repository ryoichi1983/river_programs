#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 20:20:19 2019

@author: rtsurumaki
"""

import numpy as np

class TankModelCalculator:
    """
    Ishihara & Kobatake(1979) の直列３段タンクモデルに基づいて貯留高を
    計算するクラスです。
    """
    def __init__(self, **kwargs):
        """
        タンクモデルの初期値状態設定
        ----
        Input
        ----
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

        self.alpha1 = kwargs["alpha1"]
        self.alpha2 = kwargs["alpha2"]
        self.alpha3 = kwargs["alpha3"]
        self.alpha4 = kwargs["alpha4"]
        self.beta1 = kwargs["beta1"]
        self.beta2 = kwargs["beta2"]
        self.beta3 = kwargs["beta3"]
        self.height1 = kwargs["height1"]
        self.height2 = kwargs["height2"]
        self.height3 = kwargs["height3"]
        self.height4 = kwargs["height4"]
        self.storage1 = kwargs["storage1"]
        self.storage2 = kwargs["storage2"]
        self.storage3 = kwargs["storage3"]
        self.storage = self.storage1 + self.storage2 + self.storage3

    def calculateTank(self, gridTime, rainfall):
        """
        タンク内の貯留高と流出量を計算するメソッドです。
        ----
        Input
        ----
        gridTime: 時間グリッド幅(hr)
        rainfall: 雨量(mm/hr)
        ----
        """
        """
        タンク内の貯留高と流出量を計算するメソッドです。
        ----
        Input
        ----
        gridTime: 時間グリッド幅(hr)
        rainfall: 雨量(mm/hr)
        ----
        """
        # 第 1 タンクに貯留する単位時間あたりの雨量
        self.storage1 += rainfall * gridTime

        # 第 1 タンクの単位時間あたりの流出量の定義
        if self.storage1 > self.height2:
            outflow1 = self.alpha1 * (self.storage1 - self.height1) + \
                       self.alpha2 * (self.storage1 - self.height2)
        elif self.storage1 > self.height1:
            outflow1 = self.alpha1 * (self.storage1 - self.height1)
        else:
            outflow1 = 0

        # 第 2 タンクの単位時間あたりの流出量の定義
        if self.storage2 > self.height3:
            outflow2 = self.alpha3 * (self.storage2 - self.height3)
        else:
            outflow2 = 0

        # 第 3 タンクの単位時間あたりの流出量の定義
        if self.storage3 > self.height4:
            outflow3 = self.alpha4 * (self.storage3 - self.height4)
        else:
            outflow3 = 0
            
        # 単位時間あたりの浸透量の定義
        leakage1 = self.beta1 * self.storage1
        leakage2 = self.beta2 * self.storage2
        leakage3 = self.beta3 * self.storage3

        # 第一タンクの時間ごとの貯留高の計算
        # ----------------------
        loss1 = (leakage1 + outflow1) * gridTime
        # 第一タンクの貯留高が損失高より大きい場合
        if self.storage1 >= loss1:
            self.storage1 = self.storage1 - loss1
        # 第一タンクの貯留高が損失高より小さい場合
        else:
            if self.storage1 > self.height2:
                outflow1 = self.storage1 * (self.alpha1 + self.alpha2) / \
                           (self.alpha1 + self.alpha2 + self.beta1) / gridTime
                leakage1 = self.storage1 * (self.beta1) / \
                           (self.alpha1 + self.alpha2 + self.beta1) / gridTime
            elif self.storage1 > self.height1:
                outflow1 = self.storage1 * (self.alpha1) / \
                           (self.alpha1 + self.beta1) / gridTime
                leakage1 = self.storage1 * (self.beta1) / \
                           (self.alpha1 + self.beta1) / gridTime
            else:
                outflow1 = 0 / gridTime
                leakage1 = self.storage1 / gridTime
            self.storage1 = self.storage1 - (leakage1 + outflow1) * gridTime

        # 第二タンクの時間ごとの貯留高の計算
        # ----------------------
        loss2 = (leakage2 + outflow2 - leakage1) * gridTime
        # 第二タンクの貯留高が損失高より大きい場合
        if self.storage2 >= loss2:
            self.storage2 = self.storage2 - loss2
        # 第二タンクの貯留高が損失高より小さい場合
        else:
            if self.storage2 > self.height3:
                outflow2 = self.storage2 * (self.alpha3) / \
                           (self.alpha3 + self.beta2) / gridTime
                leakage2 = self.storage2 * (self.beta2) / \
                           (self.alpha3 + self.beta2) / gridTime
            else:
                outflow2 = 0 / gridTime
                leakage2 = self.storage2 / gridTime
            self.storage2 = self.storage2 - (leakage2 + outflow2) * gridTime

        # 第三タンクの時間ごとの貯留高の計算
        # ----------------------
        loss3 = (leakage3 + outflow3 - leakage2) * gridTime
        # 第三タンクの貯留高が損失高より大きい場合
        if self.storage3 > loss3:
            self.storage3 = self.storage3 - loss3
        # 第三タンクの貯留高が損失高より少ない場合
        else:
            if self.storage3 > self.height4:
                outflow3 = self.storage3 * (self.alpha4) / \
                           (self.alpha4 + self.beta3) / gridTime
                leakage3 = self.storage3 * (self.beta3) / \
                           (self.alpha4 + self.beta3) / gridTime
            else:
                outflow3 = 0 / gridTime
                leakage3 = self.storage3 / gridTime
            self.storage3 = self.storage3 - (leakage3 + outflow3) * gridTime

        totalOutflow = outflow1 + outflow2 + outflow3
                        
        swi = self.storage1 + self.storage2 + self.storage3

        stateVector = [swi, self.storage1, self.storage2, self.storage3,
                       totalOutflow, outflow1, outflow2, outflow3,
                       leakage1, leakage2, leakage3]
        
        return stateVector

