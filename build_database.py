#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/29 10:38 下午
# @Author     : zzhen
# @File       : build_database.py
# @Software   : PyCharm
# @Description: 用于从kew.org下载数据库，并进行清洗打包，放到本地
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import argparse
import csv
import json
import os
import re
from collections import defaultdict
from urllib.error import HTTPError
from urllib.request import urlretrieve
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from Bio import SeqIO
import Easy353Lib


# detect network_connect
def network_connect() -> bool:
    connect = False
    # 网络连通 exit_code == 0
    exit_code = os.system("ping www.baidu.com")
    if not exit_code:
        connect = True
    return connect


# 解析classification.json 获取目、科、属信息
# return {"Family1":"Family","Genus1":"Genus"}
def parse_classification_json() -> dict:
    _classification_ = defaultdict(str)
    classification_file_path = os.path.join(os.path.split(
        os.path.realpath(Easy353Lib.__file__))[0], "classification.json")
    with open(classification_file_path, "r", encoding="UTF-8") as f:
        _classification_dict_ = json.load(f)
    for key, value in _classification_dict_.items():
        for i in value:
            _classification_[i] = key
    return _classification_


# 判断用户输入的拉丁名属于目、科、属的哪一级
def detect_classification(classifications: list) -> dict:
    # 设定输出结果 # {str:list[str]} eg. {"Family": ["family1", "family2"]}
    result_dict = defaultdict(list)
    # 获取拉丁学名的具体分类
    classification_dict = parse_classification_json()
    # 将输入的拉丁学名进行处理
    _classifications_ = [classification.capitalize() for classification in classifications if classification]
    for classification in _classifications_:
        tmp = classification_dict.get(classification, None)
        if tmp is None:
            print("The Latin name {} may be incorrect or the database does not have the sequences"
                  .format(classification))
        else:
            result_dict[tmp].append(classification)
    return result_dict


# 根据传入的分类信息获取需要下载的文件路径
def generate_download_info(classification_dict: dict) -> list:
    result_list = []
    kew_data_file_path = os.path.join(os.path.split(
        os.path.realpath(Easy353Lib.__file__))[0], "kew_data.csv")
    _reader_ = csv.DictReader(open(kew_data_file_path, "r", newline=""))
    # 将归属于目、科、属的数据存储在list中
    for key, value in classification_dict.items():
        # row 是存储下载信息的dict
        for row in _reader_:
            if row in result_list:
                continue
            if row[key] in value:
                result_list.append(row)
    return result_list


# 用于下载文件的函数
def download_fasta_file(_spec_info_: dict, output_dir: str):
    info = None
    url = _spec_info_['Fasta file url']
    file_path = os.path.join(output_dir, _spec_info_['Fasta file name'])
    if os.path.isfile(file_path):
        info = "INFO: File {} already exists".format(file_path)
    else:
        try:
            urlretrieve(url, file_path)
        except HTTPError:
            info = "INFO: Url {} does not exist".format(url)
        except Exception as e:
            info = "ERROR: {}".format(e)
    return info


# 根据传入的文件夹获取文件夹下所有fasta文件的路径
def generate_fasta_path(_dir_path_: str) -> list:
    # 存储所有的文件路径
    file_path_list = []
    if not os.path.isdir(_dir_path_):
        raise FileNotFoundError("The path is not a directory!")
    # 设置文件拓展名
    extension = (".fasta", ".fas", ".fa", ".faa")
    # 对文件夹进行遍历
    for _root_, _dir_, _file_ in os.walk(top=_dir_path_, topdown=False):
        if _file_:
            _file_ = [os.path.join(_root_, _one_file_) for _one_file_ in _file_ if _one_file_.endswith(extension)]
            file_path_list.extend(_file_)
    return file_path_list


# 用于根据不同物种的fasta文件，合并生成不同基因的fasta文件
def generate_gene_file(_file_path_list_: list, _output_dir_: str, _exclude_species_: list = None) -> None:
    if not os.path.isdir(_output_dir_):
        os.makedirs(_output_dir_)
    for _file_ in _file_path_list_:
        for record in SeqIO.parse(_file_, "fasta"):
            with open(os.path.join(_output_dir_, record.id + ".fasta"), "a") as _out_file_:
                try:
                    # 对于标准353基因的信息和描述的修改
                    # 使用:和空格作为分隔符
                    species_name = re.findall(".*Species:(.*)Repository.*", record.description)[0].strip().replace(" ",
                                                                                                                   "_")
                    record.id = species_name + "_" + record.id
                    if not _exclude_species_ or species_name not in _exclude_species_:
                        _out_file_.write(">" + record.id + "\n")
                        _out_file_.write(str(record.seq) + "\n")
                except IndexError:
                    print("The record {}is not standard".format(record.description))
                    _out_file_.write(">" + record.description + "\n")
                    _out_file_.write(str(record.seq) + "\n")
            _out_file_.close()
    return

def main():
    # 设置输出文件夹
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='''The script for getting the ref sequences from kew.org''')
    pars.add_argument('-i', dest="input_dir", type=str, help="The input directory of existed files from kew.org")
    pars.add_argument('-o', dest="output_dir", type=str,
                      help="The output directory stored the downloaded sequences file")
    pars.add_argument('-c', dest="classification", type=str, nargs="+",
                      help='''the classification of species that need to be downloaded''')
    pars.add_argument("-t", dest="thread", type=int, help="threads for downloading fasta file from Kew", default=4)
    pars.add_argument("-exclude", dest="exclude", type=str, nargs="+",
                      help="exclude the species that need to be downloaded")
    pars.add_argument("-exclude_file", dest="exclude_file", type=str,
                      help="The file documenting species that need to be excluded")
    pars.add_argument("-generate", dest="generate", action="store_true",
                      help="whether to generate the species that need to download")

    args = pars.parse_args()

    # 处理输入和输出文件夹
    input_dir, output_dir = args.input_dir, args.output_dir
    _thread_ = args.thread
    exclude_species = []
    if input_dir is None and output_dir is None:
        print("Input and output directory cannot both be empty!")
        exit(-1)
    if output_dir is not None:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        if input_dir is None:
            input_dir = output_dir
    if input_dir is not None and output_dir is None:
        output_dir = input_dir
    # 当需要下载文件时
    if args.classification is not None:
        print("INFO: Get species information that needs to be downloaded")
        # 根据传入的分类信息需要下载的物种信息 是元素为dict的list
        down_spec_info = generate_download_info(detect_classification(args.classification))
        # 检测已经存在的文件名
        existed_files = [os.path.basename(i) for i in generate_fasta_path(input_dir)]
        # 获取还未下载的文件
        down_spec_info = [i for i in down_spec_info if i not in existed_files]
        # 如果需要将下载的信息写入文件
        if args.generate and down_spec_info:
            print("INFO: Save the downloaded information into a csv file")
            with open(os.path.join(output_dir, "download.csv"), "a") as _out_file_:
                _writer_ = csv.DictWriter(_out_file_, down_spec_info[0].keys())
                _writer_.writeheader()
                for element in down_spec_info:
                    _writer_.writerow(element)
        # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
        # 可以使用偏函数构造出单参数函数
        print("INFO: Download species data")

        # 下载文件
        count = 1
        # 使用concurrent.features中的ThreadPoolExecutor类来实现多线程下载
        with ThreadPoolExecutor(max_workers=_thread_) as executor:
            # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
            # 可以使用偏函数构造出单参数函数 一般来说必须采用关键字参数的形式给构建出的偏函数传参，因为如果以无关键字参数的方式，该实参将试图传递给第一个形参
            # http://c.biancheng.net/view/5674.html
            futures = [executor.submit(partial(download_fasta_file, output_dir=output_dir), _spec_info_)
                       for _spec_info_ in down_spec_info]
            for future in as_completed(futures):
                if future.result() is not None:
                    print(future.result())
                if count % 10 == 0:
                    print("INFO: {} / {} has been downloaded".format(count, len(down_spec_info)))
                count += 1
    if args.exclude is not None:
        exclude_species.extend(args.exclude)
    if args.exclude_file is not None:
        if not os.path.isfile(args.exclude_file):
            print("The exclude file does not exist!")
        else:
            with open(args.exclude_file, "r") as _file_:
                _exclude_ = _file_.readlines()
                _exclude_ = [x.strip() for x in _exclude_ if x.strip()]
                exclude_species.extend(_exclude_)
    exclude_species = [species.capitalize().replace(" ", "_") for species in exclude_species]

    print("INFO: Generating species data into reference files")
    # 将input_dir下的文件和 output_dir中的文件合并
    input_file_list = generate_fasta_path(input_dir)
    output_file_list = generate_fasta_path(output_dir)
    file_path_list = list((set(input_file_list).union(set(output_file_list))))
    generate_gene_file(_file_path_list_=file_path_list,
                       _output_dir_=os.path.join(output_dir, "353gene"), _exclude_species_=exclude_species)


if __name__ == '__main__':
    main()

