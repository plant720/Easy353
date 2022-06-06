#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/2 2:21 下午
# @Author     : zzhen
# @File       : download_file.py
# @Software   : PyCharm
# @Description: 运用多进程下载353的数据,根据specimens.csv中的记录,存放在具体目,科,属的文件夹
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

from typing import Dict
from urllib.error import HTTPError
from urllib.request import urlretrieve
import os
import csv
import multiprocessing
from Bio import SeqIO
from execution_time import ExecutionTime
import re
import argparse

e = ExecutionTime(console=True)


# Generate the downloaded data table
# para: {str:list[str]} eg. {"Family": ["科1", "科2"]}
# 根据传入参数，产生需要下载的数据表
@e.timeit
def generate_data_csv(parameter: Dict, database: str = "resource/specimens.csv", data_dir: str = "data",
                      out_csv: str = "download.csv") -> None:
    with open(database, "r", newline="") as _in_file_:
        with open(os.path.join(data_dir, out_csv), "w", newline="") as _out_file_:
            _reader_ = csv.DictReader(_in_file_)
            _writer_ = csv.DictWriter(_out_file_, _reader_.fieldnames)
            _writer_.writeheader()
            key = list(parameter.keys())[0]
            # 传入的科名、属名和目名 首字母大写其他小写
            value = [i.capitalize() for i in list(parameter[key])]
            for row in _reader_:
                if row[key] in value:
                    _writer_.writerow(row)
    return


# 根据传入的文件夹获取所有fasta文件的路径
@e.timeit
def generate_fasta_path(_dir_path_: str) -> list:
    # 存储所有的文件路径
    file_path_list = []
    if not os.path.isdir(_dir_path_):
        raise FileNotFoundError("The path is not a directory")
    # 设置文件拓展名
    extension = (".fasta", ".fas", ".fa", ".fna", ".ffn", ".frn", ".faa", ".fna")
    # 对文件夹进行遍历
    for _root_, _dir_, _file_ in os.walk(top=_dir_path_, topdown=False):
        if _file_:
            _file_ = [os.path.join(_root_, _one_file_) for _one_file_ in _file_ if _one_file_.endswith(extension)]
            file_path_list.extend(_file_)
    return file_path_list


@e.timeit
# 用于根据不同物种的fasta文件，合并生成不同基因的fasta文件
def generate_gene_file(_dir_path_: str, _output_dir_: str, exclude_species: list = None) -> None:
    # 生成fasta文件的路径
    file_path_list = generate_fasta_path(_dir_path_)
    if not os.path.isdir(_output_dir_):
        os.makedirs(_output_dir_)
    for _file_ in file_path_list:
        for record in SeqIO.parse(_file_, "fasta"):
            with open(os.path.join(_output_dir_, record.id + ".fasta"), "a") as _out_file_:
                try:
                    # 对于标准353基因的信息和描述的修改
                    # 使用:和空格作为分隔符
                    species_name = re.split(r'[: ]', record.description)[4]
                    record.id = species_name + "_" + record.id
                    if exclude_species is None or species_name not in exclude_species:
                        _out_file_.write(">" + record.id + "\n")
                        _out_file_.write(str(record.seq) + "\n")
                except IndexError:
                    print("The record is not standard")
                    _out_file_.write(">" + record.description + "\n")
                    _out_file_.write(str(record.seq) + "\n")
            _out_file_.close()
    return


# 用于下载被子植物353的类
# 一般需要实例化数据csv文件和输出文件夹
class DataTool:
    # database指的是treeoflife.kew.org的数据库
    def __init__(self, data_csv: str = "resource/specimens.csv", out_dir: str = "data"):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        self.data_csv = data_csv
        self.out_dir = out_dir

    # Download the fasta from the url and save it to the data_dir
    def get_fasta_file(self, info_dict: dict, storage_by_classification: bool = False) -> None:
        fasta_dir = None
        if storage_by_classification:
            # 根据分类生成文件目录
            fasta_dir = os.path.join(self.out_dir, "fasta_file", info_dict["Order"], info_dict["Family"],
                                     info_dict["Genus"])
        else:
            fasta_dir = os.path.join(self.out_dir, "fasta_file")
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)
        # fasta_path = os.path.join(fasta_dir, ".".join(info_dict["Fasta file url"].split(".")[-3:]))
        # # 用于处理重复文件即一个物种存在两个以上的文件
        # if os.path.isfile(fasta_path):
        #     # 使用全称的名字
        #     fasta_path = os.path.join(fasta_dir, info_dict["Fasta file url"].split("/")[-1])
        fasta_path = os.path.join(fasta_dir, info_dict["Fasta file url"].split("/")[-1])
        if not os.path.isfile(fasta_path):
            try:
                urlretrieve(info_dict["Fasta file url"], fasta_path)
            except HTTPError:
                print(info_dict["Fasta file url"] + " does not exist")
        return

    # Download the data in the download_csv
    @e.timeit
    def download_data(self, storage_by_classification: bool = False, _thread_: int = 4) -> None:
        with open(self.data_csv, "r", newline="") as _in_file_:
            _reader_ = csv.DictReader(_in_file_)
            if _thread_ > 1:
                pool = multiprocessing.Pool(_thread_)
                for row in _reader_:
                    pool.apply_async(self.get_fasta_file, args=(row, storage_by_classification,))
                pool.close()
                pool.join()
            else:
                for row in _reader_:
                    self.get_fasta_file(row, storage_by_classification)
        return


# 根据传入的信息,下载文件


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Miner 
       zzhen@sculab ''')
    pars.add_argument('-csv', metavar='<str>', type=str, help='''input fasta files.''', required=True)
    pars.add_argument('-o', metavar='<str>', type=str, help='''out dir.''', required=True)
    pars.add_argument('-family', metavar='<str>', type=str, help='''the family of species that need to be downloaded''',
                      required=False, nargs="+")
    pars.add_argument('-order', metavar='<str>', type=str, help='''the order of species that need to be downloaded''',
                      required=False, nargs="+")
    pars.add_argument('-genus', metavar='<str>', type=str, help='''the genus of species that need to be downloaded''',
                      required=False, nargs="+")
    pars.add_argument("-t", metavar='<int>', type=int, help="threads for downloading fasta file from Kew", default=4)
    pars.add_argument('-s', action="store_true", help='''storage_by_classification or not''')
    pars.add_argument("-generate", action="store_true", help="generate the download csv file or not")
    pars.add_argument("-exclude", metavar='<str>', type=str, help="exclude the species that need to be downloaded",
                      required=False, nargs="+")
    args = pars.parse_args()
    download_info = {}
    if args.order is not None:
        download_info["Order"] = args.order
    if args.family is not None:
        download_info["Family"] = args.family
    if args.genus is not None:
        download_info["Genus"] = args.genus
    download_csv = args.csv
    if args.generate:
        generate_data_csv(download_info, args.csv, args.o)
        download_csv = os.path.join(args.o, "download.csv")
    data_tools = DataTool(download_csv, args.o)
    data_tools.download_data(args.s, args.t)
    generate_gene_file(args.o+"/fasta_file", args.o+"/gene_file", args.exclude)

    # python download_file.py -csv resource/specimens.csv -o ~/Desktop/test -family sapindaceae -t 4 -generate
    # para = {"Family": ["Sapindaceae"]}
    # generate_data_csv(para, "resource/specimens.csv", "/Users/zzhen/Desktop/test")
    # data_tool = DataTool("/Users/zzhen/Desktop/test/download.csv", "/Users/zzhen/Desktop/test")
    # data_tool.download_data(storage_by_classification=True, _thread_=4)
    # generate_gene_file("/Users/zzhen/Desktop/test/fasta_file/", "/Users/zzhen/Desktop/test/gene_file/")
