#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/3/28 10:22 下午
# @Author     : zzhen
# @File       : utils.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import os
import sys
import glob
import psutil
from collections import defaultdict

try:
    from reprlib import repr
except ImportError:
    pass


# def total_size(o, handlers=None, verbose=False):
#     """ Returns the approximate memory footprint an object and all of its contents.
#     Automatically finds the contents of the following builtin containers and
#     their subclasses:  tuple, list, deque, dict, set and frozenset.
#     To search other containers, add handlers to iterate over their contents:
#         handlers = {SomeContainerClass: iter,
#                     OtherContainerClass: OtherContainerClass.get_elements}
#     """
#     if handlers is None:
#         handlers = {}
#     dict_handler = lambda dst: chain.from_iterable(dst.items())
#     all_handlers = {tuple: iter,
#                     list: iter,
#                     deque: iter,
#                     dict: dict_handler,
#                     set: iter,
#                     frozenset: iter,
#                     }
#     all_handlers.update(handlers)  # user handlers take precedence
#     seen = set()  # track which object id's have already been seen
#     default_size = getsizeof(0)  # estimate sizeof object without __sizeof__
#
#     def sizeof(obj):
#         if id(obj) in seen:  # do not double count the same object
#             return 0
#         seen.add(id(obj))
#         s = getsizeof(obj, default_size)
#
#         if verbose:
#             print(s, type(obj), repr(obj), file=stderr)
#
#         for typ, handler in all_handlers.items():
#             if isinstance(obj, typ):
#                 s += sum(map(sizeof, handler(obj)))
#                 break
#         return s
#
#     return sizeof(o)


# 反向互补
def reverse_complement_all(_seq_: str) -> str:
    return _seq_.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


# 简化反向互补 速度更快
def reverse_complement_limit(_seq_: str) -> str:
    return _seq_.upper().translate(str.maketrans('ACGT', 'TGCA'))[::-1]


comp_tuple = (
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 91, 92, 93, 94, 95,
    96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127)


# 反向互补
def reverse_complement_ascii(_seq_: str) -> str:
    _seq_ = bytearray(_seq_[::-1], "ascii")
    for char in range(len(_seq_)):
        _seq_[char] = ord(comp_tuple[_seq_[char]])
    return _seq_.decode("ascii")


def log(log_file, *args):
    with open(log_file, 'a') as out:
        for i in args[:-1]:
            out.write(str(i) + ',')
        out.write(str(args[-1]) + '\n')


# 用于判断文件后缀名，用于选择具体的读取函数
def judge_type(path):
    suffix_dict = {'.gz': 0, '.fq': 1, '.fastq': 1, 'fa': 2, '.fas': 2, '.fasta': 2}
    # 设定默认值是3
    return suffix_dict.get(os.path.splitext(path)[-1].lower(), 3)


# 将字节转为字符串
# 用于处理open(file,"rb")获取的结果
def bytes_str(_input_: bytes):
    try:
        return _input_.decode('utf-8')
    except AttributeError:
        return _input_


# 如果传入的是一个文件夹，则返回该文件夹下的所有文件的文件路径列表
# 如果传入的是一个文件，则返回只包含该文件路径的列表
def get_file_list(_file_or_dir_path_: str) -> list:
    # 设置文件拓展名
    extension = (".fasta", ".fas", ".fa", ".fna", ".ffn", ".frn", ".faa", ".fna")
    _file_path_list_ = []
    if os.path.isdir(_file_or_dir_path_):
        _file_path_list_ = [_ for _ in glob.glob(os.path.join(_file_or_dir_path_, "*")) if os.path.isfile(_)
                            and _.endswith(extension)]
    elif os.path.isfile(_file_or_dir_path_):
        _file_path_list_ = [_file_or_dir_path_]
    else:
        pass
    return _file_path_list_


# 返回文件大小 以MB为单位
def get_file_size(file_path: str) -> float:
    fsize = os.path.getsize(file_path)
    fsize = fsize / float(1024 * 1024)
    return round(fsize, 2)


# 判断文件中的reads数量及read长度
def get_reads_info(file_path: str) -> tuple:
    reads_num, reads_len = 0, 0
    # 读取fasta数据
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                reads_num += 1
            else:
                reads_len += len(line.strip().upper().replace("N", ""))
    return reads_num, reads_len


def get_seq_avg_length(file_path: str) -> float:
    reads_num, reads_len = get_reads_info(file_path)
    return reads_len / reads_num


# 用于返回与参考序列相关的信息:
# 基因名/文件名 与其对应的文件路径
# 基因名/文件名 与其对应的参考序列的平均长度
# todo:这一步实际上已经读取了一边参考序列,为了减少操作,可以在make_ref_kmer_dict这一个函数中返回
def get_ref_info(reference_path: str) -> tuple:
    # 根据reference路径获取_ref_path_list_
    _ref_path_list_ = get_file_list(reference_path)
    if not _ref_path_list_:
        print("Error: reference file is not exist, please check your input.")
        sys.exit(-1)
    # 记录reference文件的read平均长度
    _ref_length_dict_ = defaultdict(int)
    # 记录reference文件的路径信息
    _ref_path_dict_ = defaultdict(str)
    for ref_file_path in _ref_path_list_:
        if judge_type(ref_file_path) == 2:
            # 通过参考文件的文件名来标定不同的参考基因
            file_name_gene = os.path.basename(ref_file_path).split('.')[0]
            # 参考序列的路径字典
            _ref_path_dict_[file_name_gene] = ref_file_path
            ref_seq_count, _ref_length_dict_[file_name_gene] = get_reads_info(ref_file_path)
            # 计算平均长度
            _ref_length_dict_[file_name_gene] = int(_ref_length_dict_[file_name_gene] / max(ref_seq_count, 1))
    return _ref_length_dict_, _ref_path_dict_


# 创建哈西字典,将reference存入字典中
# _ref_reverse_complement_会影响建立的哈西字典的大小
# 如果要节约内存，将_ref_reverse_complement_设置为False,那么在过滤测序数据时，需要对read进行reverse complement
def make_ref_kmer_dict(reference_path: str, _kmer_size_: int, _ref_reverse_complement_: bool = False,
                       _pos_: bool = True, _print_: bool = True, ref_number: int = None) -> dict:
    # 从reference文件中获取文件路径列表
    files_list = get_file_list(reference_path)
    if not files_list:
        print('Reference is invalid. Please check the path!')
        sys.exit(1)
    ref_kmer_dict = defaultdict(list)
    gene_number_count = 0
    for file in files_list:
        # 用于记录每个基因的序列数量
        ref_count = 0
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        file_name_gene = os.path.basename(file).split(".")[0]
        # 跳过第一行>
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            gene_number_count += 1
            # 设定单个参考基因的最大序列数量
            ref_count += 1
            if ref_number is not None:
                if ref_count >= ref_number:
                    break
            # 保留seq序列中的数字和字母 去除换行符
            ref_seq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            for j in range(0, len(ref_seq) - _kmer_size_ + 1):
                kmer_info_list, ref_kmer = [], ref_seq[j:j + _kmer_size_]
                if ref_kmer in ref_kmer_dict:
                    kmer_info_list = ref_kmer_dict[ref_kmer]
                    # temp_list[0]是ref_kmer出现次数
                    kmer_info_list[0] += 1
                    # 如果需要添加位置信息
                    if _pos_:
                        # temp_list[1]是ref_kmer出现的百分比位置总和
                        kmer_info_list[1] += (j + 1) / len(ref_seq)
                    if file_name_gene not in kmer_info_list:
                        # temp_list[2]是ref_kmer出现的参考基因基因列表
                        kmer_info_list.append(file_name_gene)
                else:
                    kmer_info_list = [1, (j + 1) / len(ref_seq), file_name_gene] if _pos_ else [1, 0,
                                                                                                file_name_gene]
                    # kmers字典的key是kmer_seq value是一个字典, list[0] 是kmer出现次数
                    # list[1] 是kmer出现位置 list[2:]是该kmer在哪些文件(基因)中出现
                    # sys.intern(ref_kmer)用于加快内存
                    ref_kmer_dict[sys.intern(ref_kmer)] = kmer_info_list
            # _ref_reverse_complement_ 对reference中的序列进行反向互补
            if _ref_reverse_complement_:
                ref_seq = reverse_complement_limit(ref_seq)
                for j in range(0, len(ref_seq) - _kmer_size_ + 1):
                    kmer_info_list, ref_kmer = [], ref_seq[j:j + _kmer_size_]
                    if ref_kmer in ref_kmer_dict:
                        kmer_info_list = ref_kmer_dict[ref_kmer]
                        kmer_info_list[0] += 1
                        if _pos_:
                            # 反向序列的pos值为负值
                            kmer_info_list[1] -= (j + 1) / len(ref_seq)
                        if file_name_gene not in kmer_info_list:
                            kmer_info_list.append(file_name_gene)
                    else:
                        kmer_info_list = [1, -(j + 1) / len(ref_seq), file_name_gene] if _pos_ else [1, 0,
                                                                                                     file_name_gene]
                        ref_kmer_dict[sys.intern(ref_kmer)] = kmer_info_list
            # json.dumps(my_dictionary)后面进行检测
            # print(json.dumps(ref_kmer_dict))
            if _print_:
                print('Mem.:', round(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, 2),
                      'G, Num. of Seq:', gene_number_count, end='\r')
            else:
                if gene_number_count % 1000 == 0:
                    print("INFO: Mem. used: {} G, Num. of Seq: {}".format(round(
                        psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, 2), gene_number_count))
        infile.close()
    return ref_kmer_dict


def dynamic_limit_v2(_read_kmer_dict_: dict, smoother_level: int = 4, smoother_size: int = 64):
    count_list = [0] * smoother_size
    for x in _read_kmer_dict_:
        if _read_kmer_dict_[x][0] <= smoother_size:
            count_list[_read_kmer_dict_[x][0] - 1] += 1
    # 平滑器
    for i in range(smoother_level):
        for x in range(1, smoother_size - smoother_level + i):
            if count_list[x + smoother_level - 1 - i] < count_list[x - 1] \
                    and count_list[x] < count_list[x + smoother_level - i]:
                return x + 1
    return 2


# 动态limit
def dynamic_limit(_read_kmer_dict_: dict, gen_avg_len: int, ref_rate: int = 1.5, list_size: int = 256) -> int:
    count_list = [0] * list_size
    for _read_kmer_ in _read_kmer_dict_:
        # _read_kmer_dict_[x][0]是kmer出现次数 即计算出现出现n次的所有kmer的总数
        if _read_kmer_dict_[_read_kmer_][0] <= list_size:
            count_list[_read_kmer_dict_[_read_kmer_][0] - 1] += 1
    # 计算器
    F0, sum_f = len(_read_kmer_dict_), 0
    for x in range(list_size):
        sum_f += count_list[x]
        if (F0 - sum_f) / 2 < gen_avg_len * ref_rate:
            return max(2, x)
    return 2
