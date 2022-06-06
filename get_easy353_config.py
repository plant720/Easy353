#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/29 8:50 下午
# @Author     : zzhen
# @File       : get_easy353_config.py
# @Software   : PyCharm
# @Description:
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import platform
import sys

# system info
SYSTEM_NAME = ""
if platform.system() == "Linux":
    SYSTEM_NAME = "linux"
elif platform.system() == "Darwin":
    SYSTEM_NAME = "macOS"
elif platform.system() == "Windows":
    SYSTEM_NAME = "windows"
else:
    sys.stdout.write("Error: currently Easy353 is not supported for " + platform.system() + "! ")
    sys.exit()

# python version
MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    pass
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    pass
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)






