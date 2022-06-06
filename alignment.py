#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:59 下午
# @Author     : zzhen
# @File       : alignment.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import gc
import platform
import os
import glob
import shutil


# Judge the system and select the muscle version
def get_muscle():
    if platform.system() == 'Linux':
        return 'lib/muscle/muscle3.8.31_i86linux64'
    elif platform.system() == 'Windows':
        return 'lib/muscle/muscle3.8.31_i86win32.exe'
    elif platform.system() == 'Darwin':
        return 'lib/muscle/muscle3.8.31_i86darwin64'
    else:
        print('Error: Unsupported system!')
        return None


# 返回文件夹下的所有基因文件
def get_file_dict(_dir_path_: str) -> dict:
    _file_path_dict_ = {}
    if os.path.isdir(_dir_path_):
        for _ in glob.glob(os.path.join(_dir_path_, '*')):
            ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
            if os.path.isfile(_) and _.endswith(ext):
                _file_path_dict_[_.split('/')[-1].split(".")[0]] = _
    else:
        print('Error: The input path is not a directory!')
    return _file_path_dict_


if __name__ == '__main__':
    # 初始化参数
    # 设置输出文件夹
    _out_dir_ = "/Users/zzhen/Desktop/test"
    reference_dir = "/Users/zzhen/Desktop/assembly/test-data/353_order/Apiales"
    # 是否对每个拼接结果进行alignment 默认是只对assemble成功的结果进行alignment
    _all_ = False

    combine_dir = os.path.join(_out_dir_, "combine")
    alignment_dir = os.path.join(_out_dir_, "alignment")
    if not os.path.isdir(combine_dir):
        os.mkdir(combine_dir)
    if not os.path.isdir(alignment_dir):
        os.mkdir(alignment_dir)
    # 获取参考文件路径
    reference_file_dict = get_file_dict(reference_dir)
    assemble_file_dict = get_file_dict(os.path.join(_out_dir_, 'contig'))
    if _all_:
        assemble_file_dict.update(get_file_dict(os.path.join(_out_dir_, 'short_contig')))

    for gen_name, gene_file in assemble_file_dict.items():
        print("combine:", gen_name, gene_file)
        if gen_name in reference_file_dict.keys():
            shutil.copyfile(reference_file_dict[gen_name], os.path.join(combine_dir, gen_name + ".combine.fasta"))
            with open(os.path.join(combine_dir, gen_name + ".combine.fasta"), 'a') as f:
                f.writelines(open(gene_file).readlines())

    print("combine done")

    combine_dir_file_dict = get_file_dict(combine_dir)
    for gen_name, gene_file in combine_dir_file_dict.items():

        print("alignment:", gen_name, gene_file)
        # muscle -in test.fasta -out test.alignment.fasta
        os.system("{} -in {} -out {}".format(get_muscle(), gene_file,
                                             os.path.join(alignment_dir, gen_name + ".alignment.fasta")))
        print("alignment done")
