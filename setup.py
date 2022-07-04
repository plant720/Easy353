#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/6/7 10:39 上午
# @Author     : zzhen
# @File       : setup.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

"""Setup the Easy353 environment."""

# Test pip
# 1) Clean the /dist directory
# 2) python3 setup.py sdist bdist_wheel
# 3) pip install --index-url https://test.pypi.org/simple/
#    --extra-index-url https://pypi.org/simple atram
# 4) twine upload --repository-url https://test.pypi.org/legacy/ dist/*


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

setup(
    name="Easy353",
    version="1.5.0",
    author="zzhen",
    author_email="zzhen0302@163.com",
    description="""a tool for recovering the Angiosperms353 gene set(AGS), 
                   which can filter and de novo assemble reads from sequencing data,
                   helping users capture AGS accurately and effectively.""",
    license="MIT License",
    # 项目主页
    url="https://github.com/plant720/Easy353",
    python_requires='>=3.7',
    install_requires=install_dependencies,
    entry_points={
        'console_scripts': [
            'easy353 = easy353:main',
            'build_reference = build_database:main',
        ],
    },
    packages=["Easy353Lib"],
    package_data={'Easy353Lib': ["Easy353Lib/kew_data.csv", "Easy353Lib/classification.json", "Easy353Lib/genes.csv"]},
    include_package_data=True,
    zip_safe=False
)
