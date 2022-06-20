#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/6/13 5:21 下午
# @Author     : zzhen
# @File       : parse_kew.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import requests
from bs4 import BeautifulSoup
import collections
import csv
import json

# kew fasta文件存储位置
URL = "http://sftp.kew.org/pub/paftol/current_release/fasta/by_recovery/"


def parse_html() -> dict:
    resource = requests.get(URL)
    soup = BeautifulSoup(resource.text, "html.parser")
    species_info = collections.defaultdict(dict)
    # 获取网页上的文本的url
    for element in soup.find_all('a', href=True):
        fasta_file_url = element.get('href')
        fasta_file_name = element.get_text()
        if fasta_file_name.endswith(".fasta"):
            species_info[fasta_file_name] = URL + fasta_file_url
    return species_info


if __name__ == "__main__":
    parse_html()




