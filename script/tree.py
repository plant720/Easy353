#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 3:00 下午
# @Author     : zzhen
# @File       : tree.py
# @Software   : PyCharm
# @Description: 根据给出的多序列比对文件,生成tree文件
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.


# here put the import lib
import os
import subprocess
import sys
import glob
import platform


# Judge the system and select the RaxML version
# todo:选择或者编译指定版本的RaxML
def get_RaxML():
    if platform.system() == 'Linux':
        return 'lib/RaxML/raxmlHPC-SSE3'
    elif platform.system() == 'Darwin':
        return 'lib/RaxML/raxmlHPC-SSE3'
    elif platform.system() == 'Windows':
        return 'lib/RaxML/raxmlHPC-SSE3'
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


def make_tree(_alignment_dir_: str, _out_dir_: str, _out_group_: str, getname):
    alignment_file_dict = get_file_dict(_alignment_dir_)
    # 通过RaxML进行建树
    alignment_count = 0  # 记录比对文件个数
    for gen_name, alignment_file in alignment_file_dict.items():
        alignment_count += 1
        print("Going to make_tree_withRaxML: {}      [Schedule: {}/{}]".format(alignment_file,
                                                                               alignment_count,
                                                                               len(alignment_file_dict)))
        print("---------------------")

        # 创建树文件夹
        input_file = '''"{}"'''.format(alignment_file)
        output_dir = os.path.join(_out_dir_, gen_name)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_dir = '''"{}"'''.format(output_dir)

        # 创建树文件
        cmd = "{} -# 10000 -p 12345 -m GTRGAMMA -n \
        tre -s {} -o {} -w {}".format(get_RaxML(), input_file, _out_group_, output_dir)
        subprocess.call(cmd, shell=True)

    print("*** [make_tree_withRaxML] Done.")


if __name__ == '__main__':
    dir_name = sys.argv[1]
    type = sys.argv[2]
    outgroup = sys.argv[3]
    getname = sys.argv[4]
    make_tree(dir_name, type, outgroup, getname)
