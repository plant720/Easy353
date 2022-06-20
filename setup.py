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

import re
from setuptools import setup, find_packages


def readme():
    """Get README.md content."""
    with open("README.md", 'r') as f:
        return f.read()


def license_():
    """Get LICENSE.txt content."""
    with open("LICENSE", 'r') as f:
        return f.read()


def find_version():
    """Read version from db.py."""
    regex = r"^ATRAM_VERSION = ['\"]v?([^'\"]*)['\"]"
    with open("./lib/db.py", 'r') as f:
        match = re.search(regex, f.read(), re.M)
        if match:
            return match.group(1)

    raise RuntimeError("Unable to find version string.")


def find_requirements():
    """Read requirements.txt file and returns list of requirements."""
    with open("requirements.txt", 'r') as f:
        return f.read().splitlines()


setup(
    name="Easy353",
    version=find_version(),
    author="zzhen",
    author_email="zzhen0302@gmail.com",
    packages=find_packages(),
    install_requires=find_requirements(),
    description="""a tool for recovering the Angiosperms353 gene set(AGS), 
                   which can filter and de novo assemble reads from sequencing data,
                   helping users capture AGS accurately and effectively.""",
    long_description=readme(),
    license=license_(),
    # 项目主页
    url="https://github.com/plant720/Easy353",
    python_requires='>=3.6',
    scripts=[]
)
