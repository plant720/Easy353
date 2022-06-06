#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/3 9:28 下午
# @Author     : zzhen
# @File       : process_data.py
# @Software   : PyCharm
# @Description: 将下载的353数据进行处理,按照基因名进行存储
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import os
from Bio import Seq, SeqIO
import re


def dir_walk(_dir_path_: str) -> list:
    file_list = []
    if not os.path.isdir(_dir_path_) and os.path.isfile(_dir_path_):
        file_list.append(_dir_path_)
    for _root_, _dir_, _file_ in os.walk(top=_dir_path_, topdown=False):
        if _file_:
            _file_ = [os.path.join(_root_, _one_file_) for _one_file_ in _file_]
            file_list.extend(_file_)
    return file_list


# 传入一个文件夹，根据文件夹将单个物种的353基因存放在对应的基因文件中
def split_353_genus_files(_dir_path_: str, _out_dir_root_path_: str):
    _out_dir_path_ = os.path.join(_out_dir_root_path_, _dir_path_)
    if not os.path.exists(_out_dir_path_):
        os.makedirs(_out_dir_path_)
    _353_file_list = [os.path.join(_dir_path_, i) for i in os.listdir(_dir_path_) if i != ".DS_Store"]
    for _353_file in _353_file_list:
        for record in SeqIO.parse(_353_file, "fasta"):
            gene_code = record.id
            try:
                species_name = re.split(r'[: ]', record.description)[4]
                record.id = species_name + "_" + record.id
                record.name = species_name + "_" + record.name
                record.description = species_name + "_" + record.description
                # 修改SeqIO的write函数,将w修改为a+
                SeqIO.write(record, os.path.join(_out_dir_path_, gene_code + ".fasta"), "fasta")
            except IndexError:
                print(_353_file, ":", record.description)


# 传入一个文件夹，根据文件夹将单个物种的353基因存放在对应的基因文件中
def split_353_family_files(_dir_path_: str, _out_dir_root_path_: str):
    _out_dir_path_ = os.path.join(_out_dir_root_path_, _dir_path_)
    if not os.path.exists(_out_dir_path_):
        os.makedirs(_out_dir_path_)
    _353_file_list = []
    for _root_, _, _file_ in os.walk(_dir_path_):
        if _file_:
            _353_file_list.extend([os.path.join(_root_, i) for i in _file_ if i != ".DS_Store"])
    # print(_353_file_list)
    for _353_file in _353_file_list:
        for record in SeqIO.parse(_353_file, "fasta"):
            # 修改SeqIO的write函数,将w修改为a+
            SeqIO.write(record, os.path.join(_out_dir_path_, os.path.basename(_353_file)), "fasta")


if __name__ == '__main__':
    # 用来处理关于属的数据
    # print("------开始进行属的处理------")
    # genus_dir = []
    # for _root_, _, _ in os.walk(top="/Users/zzhen/Desktop/assembly/test-data/353", topdown=False):
    #     if len(_root_.split("/")) == 10:
    #         genus_dir.append(_root_.replace("/Users/zzhen/Desktop/assembly/test-data/353/", ""))
    # print(genus_dir)
    #
    # # _dir_path_ = "Acorales/Acoraceae/Acorus"
    # # _out_dir_root_path_ = "test"
    # # _out_dir_path_ = os.path.join(_out_dir_root_path_, _dir_path_)
    # # if not os.path.exists(_out_dir_path_):
    # #     os.makedirs(_out_dir_path_)
    # # _353_file_list = [os.path.join(_dir_path_, i) for i in os.listdir(_dir_path_)]
    # # for _353_file in _353_file_list:
    # #     for record in SeqIO.parse(_353_file, "fasta"):
    # #         gene_code = record.id
    # #         species_name = re.split(r'[: ]', record.description)[4]
    # #         record.id = species_name + "_" + record.id
    # #         record.name = species_name + "_" + record.name
    # #         record.description = species_name + "_" + record.description
    # #         # 修改SeqIO的write函数,将w修改为a+
    # #         SeqIO.write(record, os.path.join(_out_dir_path_, gene_code + ".fasta"), "fasta")

    # _project_dir_ = "/Users/zzhen/Desktop/assembly/test-data/353"
    # os.chdir(_project_dir_)
    # count = 0
    # for _genus_dir_path_ in genus_dir:
    #     split_353_genus_files(_genus_dir_path_, "result")
    #     count += 1
    #     if count % 100 == 0:
    #         print("已进行{}/{}".format(count, len(genus_dir)))
    # print("------属已处理结束------")

    #   用于处理科的数据
    # print("------开始进行科的处理------")
    # family_dir = []
    # for _root_, _, _ in os.walk(top="/Users/zzhen/Desktop/assembly/test-data/353_genus", topdown=False):
    #     if len(_root_.split("/")) == 9:
    #         family_dir.append(_root_.replace("/Users/zzhen/Desktop/assembly/test-data/353_genus/", ""))
    # print(family_dir)
    #
    # _project_dir_ = "/Users/zzhen/Desktop/assembly/test-data/353_genus"
    # os.chdir(_project_dir_)

    # # 传入科的路径
    # _dir_path_ = "Alismatales/Zosteraceae"
    # _out_dir_root_path_ = "family_results"
    # _out_dir_path_ = os.path.join(_out_dir_root_path_, _dir_path_)
    # if not os.path.exists(_out_dir_path_):
    #     os.makedirs(_out_dir_path_)
    # _353_file_list = []
    # for _root_, _, _file_ in os.walk(_dir_path_):
    #     if _file_:
    #         _353_file_list.extend([os.path.join(_root_, i) for i in _file_ if i != ".DS_Store"])
    # # print(_353_file_list)
    # for _353_file in _353_file_list:
    #     for record in SeqIO.parse(_353_file, "fasta"):
    #         # 修改SeqIO的write函数,将w修改为a+
    #         SeqIO.write(record, os.path.join(_out_dir_path_, os.path.basename(_353_file)), "fasta")

    # count = 0
    # for _family_dir_path_ in family_dir:
    #     split_353_family_files(_family_dir_path_, "family_result")
    #     count += 1
    #     if count % 100 == 0:
    #         print("已进行{}/{}".format(count, len(family_dir)))
    # print("------科已处理结束------")

    # 用来处理目的数据
    print("------开始进行目的处理------")
    order_dir = []
    for _root_, _, _ in os.walk(top="/Users/zzhen/Desktop/assembly/test-data/353_family", topdown=False):
        if len(_root_.split("/")) == 8:
            order_dir.append(_root_.replace("/Users/zzhen/Desktop/assembly/test-data/353_family/", ""))
    print(order_dir)

    _project_dir_ = "/Users/zzhen/Desktop/assembly/test-data/353_family"
    os.chdir(_project_dir_)

    count = 0
    for _order_dir_path_ in order_dir:
        split_353_family_files(_order_dir_path_, "order_result")
        count += 1
        print("已进行{}/{}".format(count, len(order_dir)))
    print("------目已处理结束------")

