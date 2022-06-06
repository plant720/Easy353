#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/1 10:27 下午
# @Author     : zzhen
# @File       : ref_dir_process.py
# @Software   : PyCharm
# @Description: 将传入文件处理为hybpiper可以读取的格式
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

from Bio import SeqIO
import os
import re


def deal_ref_for_hybpiper(_input_dir_: str, _output_file_: str):
    _file_list_ = [os.path.join(_input_dir_, _file_) for _file_ in os.listdir(_input_dir_) if _file_.endswith('.fasta')]
    _out_file_ = open(_output_file_, 'a')
    a = []
    for _file_ in _file_list_:
        for record in SeqIO.parse(_file_, "fasta"):
            species_name = record.id.rsplit("_", 1)[0]
            gene_id = record.id.rsplit("_", 1)[1]
            if species_name+gene_id in a:
                species_name = species_name + "_2"
            # record_info = record.id.split("_", 1)
            # info = record_info[1] + "-" + record_info[0]
            _out_file_.write(">" + species_name + "-" + gene_id + "\n")
            _out_file_.write(str(record.seq) + "\n")
            a.append(species_name+gene_id)
    _out_file_.close()


# 根据传入的文件夹和需要去除的物种列表，处理文件
def deal_ref_for_easy353(_input_dir_: str, exclude_species: list, _output_dir_: str):
    _file_list_ = [os.path.join(_input_dir_, _file_) for _file_ in os.listdir(_input_dir_) if _file_.endswith('.fasta')]
    if not os.path.isdir(_output_dir_):
        os.makedirs(_output_dir_)
    for _file_ in _file_list_:
        _file_name_ = os.path.basename(_file_)
        _out_file_ = open(os.path.join(_output_dir_, _file_name_), 'a')
        for record in SeqIO.parse(_file_, "fasta"):
            record_info = record.id.split("_")
            species_name = " ".join(record_info[:-1])
            if species_name not in exclude_species:
                _out_file_.write(">" + record.id + "\n")
                _out_file_.write(str(record.seq) + "\n")
        _out_file_.close()


# 将一个物种的353基因序列提取出来，按照基因名存储
def split_species_file_to_gene_dir(_input_file_: str):
    _output_dir_ = _input_file_.replace(".a353.fasta", "_standard")
    if not os.path.isdir(_output_dir_):
        os.makedirs(_output_dir_)
    for record in SeqIO.parse(_input_file_, "fasta"):
        gene_code = record.id
        species_name = re.split(r'[: ]', record.description)[4]
        with open(os.path.join(_output_dir_, gene_code + ".fasta"), 'a') as _out_file_:
            _out_file_.write(">" + species_name + "_" + record.id + "\n")
            _out_file_.write(str(record.seq) + "\n")


if __name__ == '__main__':
    # species_name = ['Sinocarum schizopetalum var. bijiangense', 'Sinocarum cruciatum', 'Sinocarum cruciatum var.
    # linearilobum', 'Sinocarum filicinum', 'Sinocarum vaginatum', 'Sinocarum schizopetalum', 'Sinocarum
    # schizopetalum', 'Acronema sichuanense', 'Acronema schneideri']
    # deal_ref_for_easy353("/Users/zzhen/Desktop/Apiaceae", species_name, "/Users/zzhen/Desktop/Apiaceae_easy353")

    # species_name = ["Aegopodium podagraria", "Anthriscus sylvestris", "Chamaesium paradoxum",
    #                 "Cyclospermum leptophyllum", "Daucus carota", "Foeniculum vulgare",
    #                 "Haplosphaera phaea", "Oenanthe javanica", "Ostericum grosseserratum"]
    # deal_ref_for_easy353("/Users/zzhen/Desktop/Apiaceae", species_name, "/Users/zzhen/Desktop/Apiaceae_easy353")
    # filelist = [os.path.join("/Users/zzhen/Desktop/ref", _file_) for _file_ in os.listdir("/Users/zzhen/Desktop/ref") if
    #             _file_.endswith('.fasta')]
    # for filelist in filelist:
    #     split_species_file_to_gene_dir(filelist)
    deal_ref_for_hybpiper("/Users/zzhen/Desktop/ref_Apiaceae_9", "/Users/zzhen/Desktop/Apiaceae_hybpiper.fasta")
