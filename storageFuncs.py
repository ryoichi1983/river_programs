# -*- coding: utf-8 -*-
"""
Created on Tue Jul 03 15:12:00 2019

@author: ryoichi.tsurumaki

"""

import math
import numpy as np


def classicOneValueStorageFunc(t, y, a, p, rainfall):
    """
    古典的一価非線型貯留関数 (有効雨量を考慮)
    s = a q^p
    ds / dt = r - q - bq
    -> dy / dt = (r - (1 + b) * y^1/p) / a, where y = q^p
    貯留量, q: 流量, t: 時間, r: 雨量
    parameters: p, a, b
    """
    dydt = (rainfall - y**(1/p)) / a  # 浸透なし
    return dydt


def classicTwoValueStorageFunc(t, y, a, b, m, n, rainfall):
    """
    古典的二価非線型貯留関数 (有効雨量を考慮)
    s = a * q^m +  b * d(q^n) / dt
    ds / dt = r - q
    -> (1) dy/dt = -(am/bn) * x^(m/n - 1) * y - (1 + c) * x^(1/n) / b + r / b
       (2) dx/dt = y, where x = q^n, 
    s: 貯留量, q: 流量, t: 時間, r: 雨量
    parameters: a, b, m, n
    """
    dxdt = y[1]
    dydt = - a/b * m/n * y[0]**(m/n - 1) * y[1] - y[0]**(1/n) / b + rainfall / b  # 浸透なし
    return [dxdt, dydt]  # (x, y)


def twoStepTwoValueStorageFunc(t, y, k11, k12, k13, k21, k22, p1, p2, rainfall):
    """
    二段タンク型二価非線型貯留関数
        * タンク1：二価非線型貯留関数
        * タンク2：二価線型貯留関数
    s1 = k11 * q1^p1 + k12 * d(q1^p2)/dt
    ds1/dt = r - q1 - k13 * q1
    s2 = k21 * q2 +  k22 * d(q2)/dt
    ds2/dt = k13 * q1 - q2
    -> (1) dy1/dt = y2, where y1 = q1^p2
       (2) dy2/dt = - (k11/k12) * (p1/p2) * y1^(p1/p2 - 1) * y2 - (1 + k13) * y1^(1/p2) / k12 + r / k12
       (3) dy3/dt = y4, where y3 = q2
       (4) dy4/dt = - (k21/k22) * y4 - y3 / k22 + k13 * y1^(1/p2) / k22
    si: タンク i の貯留量, qi: タンク i の流量, t: 時間, r: 雨量
    parameters: k11, k12, k13, k21, k22, p1, p2
    """
    dy1dt = y[1]
    dy2dt = - k11/k12 * p1/p2 * y[0]**(p1/p2-1) * y[1] - \
              (1 + k13)/k12 * y[0]**(1/p2) + rainfall/k12  # 浸透あり
    dy3dt = y[3]
    dy4dt = - k21/k22 * y[3] - y[2]/k22 + k13/k22 * y[0]**(1/p2)  # 浸透なし
    return [dy1dt, dy2dt, dy3dt, dy4dt]  # (y1, y2, y3, y4)
