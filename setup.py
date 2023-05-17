#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/5/16 10:39 上午
# @Author     : zzhen
# @File       : setup.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2023 by sculab, All Rights Reserved.
import platform
import sys
import setuptools
from setuptools import setup

# python libs
install_dependencies = []
try:
    import biopython
except ImportError:
    install_dependencies.append("biopython>=1.70")
else:
    sys.stdout.write("Existed module numpy " + str(biopython.__version__) + "\n")

try:
    import psutil
except ImportError:
    install_dependencies.append("psutil>=5.6.6")
else:
    sys.stdout.write("Existed module sympy " + str(psutil.__version__) + "\n")
try:
    import requests
except ImportError:
    install_dependencies.append("requests[security]")
else:
    sys.stdout.write("Existed module requests " + str(requests.__version__) + "\n")
try:
    from bs4 import beautifulsoup4
except ImportError:
    install_dependencies.append("beautifulsoup4>=4.9.0")
else:
    sys.stdout.write("Existed module beautifulsoup4 " + str(beautifulsoup4.__version__) + "\n")

sys.stdout.write("Python " + str(sys.version).replace("\n", " ") + "\n")
sys.stdout.write("PLATFORM: " + " ".join(platform.uname()) + "\n")
sys.stdout.write("Using setuptools " + str(setuptools.__version__) + "\n")

scripts_to_install = ["easy353.py", "build_database.py"]

setup(
    name="Easy353",
    version="2.0.1",
    author="zzhen",
    author_email="zzhen0302@163.com",
    description="""a tool specifically designed to recover the Angiosperms353 gene set (AGS). It effectively filters AGS-related reads from high-throughput sequencing data, and accurately recovers AGS using its optimized reference-guided assembler.""",
    license="MIT License",
    url="https://github.com/plant720/Easy353",
    python_requires='>=3.7',
    install_requires=install_dependencies,
    scripts=scripts_to_install,
    packages=["Easy353Lib"],
    package_data={'Easy353Lib': ["Easy353Lib/kew_data.csv", "Easy353Lib/classification.json", "Easy353Lib/genes.csv"]},
    include_package_data=True,
    zip_safe=False
)
