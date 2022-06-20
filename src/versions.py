#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/6/7 2:59 下午
# @Author     : zzhen
# @File       : versions.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

def get_version():
    return versions[0]["number"]


versions = [
    {
        "number": "1.2.0.1",
        "features": ["1: better log",
                     "2: added new strategy for picking seeds ",
                     "3: fix bug of error assembly"],
        "date": "2022-06-07",
    }
]
