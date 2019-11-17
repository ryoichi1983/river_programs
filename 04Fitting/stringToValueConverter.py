#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 21:11:33 2019

@author: rtsurumaki
"""

import numpy as np

#--------------------------------------------------------------------------
# #  データの型を文字列から数値に変換するクラスの定義
#--------------------------------------------------------------------------
class StringToValueConverter:
    def convert(self, item_str):

        # . と - と半角空白を除く文字列が数字ならば float 型に変換します。
        if item_str.replace('.','').replace('-','').replace(' ','').isdigit():
            return float(item_str)
    
        # 空のデータならば空の文字列で返します。
        if not item_str:
            return ''
    
        # データなし（///）の場合 NaN を返す。
        if item_str == '///':
            return np.NaN

        # 現象なし（-）の場合 NaN を返す。
        elif item_str == "-":
            return np.NaN
        
        return item_str
