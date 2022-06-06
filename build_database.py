#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/29 10:38 下午
# @Author     : zzhen
# @File       : build_database.py
# @Software   : PyCharm
# @Description: 用于从kew.org下载数据库，并进行清洗打包，放到本地
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
from urllib.error import HTTPError
from urllib.request import urlretrieve
import requests
from bs4 import BeautifulSoup
import collections
import os
import multiprocessing


def get_fasta_file(url, file_path) -> None:
    if os.path.isfile(file_path):
        return
    try:
        urlretrieve(url, file_path)
    except HTTPError:
        print(url + " does not exist")
    return


if __name__ == '__main__':
    # 设置输出文件夹
    out_dir = "/Users/zzhen/Desktop/easy353"
    URL = "http://sftp.kew.org/pub/paftol/current_release/fasta/by_recovery/"
    # 设置进程数
    processes = 4
    resource = requests.get(URL)
    soup = BeautifulSoup(resource.text, "html.parser")
    species_info = collections.defaultdict(dict)
    for element in soup.find_all('a', href=True):
        link = element.get('href')
        txt = element.get_text()
        if txt.endswith(".fasta"):
            species_info[txt] = URL + link

    print(len(species_info))
    pool = multiprocessing.Pool(processes=processes)
    for key, value in species_info.items():
        pool.apply_async(get_fasta_file, args=(value, os.path.join(out_dir, key),))
    pool.close()
    pool.join()
