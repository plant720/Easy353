#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/3/28 10:22 下午
# @Author     : zzhen
# @File       : utils.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import os
import random
import sys
import glob
import psutil
from collections import defaultdict


# reverse and complement
def reverse_complement_all(_seq_: str) -> str:
    return _seq_.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


# Simplifying reverse complement, without considering the concatenated bases
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


# Used to determine the file suffix name, used to select the specific read function
def judge_type(path):
    suffix_dict = {'.gz': 0, '.fq': 1, '.fastq': 1, 'fa': 2, '.fas': 2, '.fasta': 2}
    # 设定默认值是3
    return suffix_dict.get(os.path.splitext(path)[-1].lower(), 3)


# Converting bytes to strings
# Used to process the result obtained by open(file, "rb")
def bytes_str(_input_: bytes):
    try:
        return _input_.decode('utf-8')
    except AttributeError:
        return _input_


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


# get file size: MB
def get_file_size(file_path: str) -> float:
    fsize = os.path.getsize(file_path)
    fsize = fsize / float(1024 * 1024)
    return round(fsize, 2)


#  Determine the number of reads in the file and the length of the read
def get_reads_info(file_path: str) -> tuple:
    reads_num, reads_len = 0, 0
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


# Used to return information related to the reference sequence:
# Gene name/filename and its corresponding file path
# The average length of the gene/file name and its corresponding reference sequence
def get_ref_info(reference_path: str) -> tuple:
    _ref_path_list_ = get_file_list(reference_path)
    if not _ref_path_list_:
        print("Error: reference file is not exist, please check your input.")
        sys.exit(-1)
    _ref_length_dict_ = defaultdict(int)
    _ref_path_dict_ = defaultdict(str)
    for ref_file_path in _ref_path_list_:
        if judge_type(ref_file_path) == 2:
            # identify different reference genes by the file name of the reference file
            file_name_gene = os.path.basename(ref_file_path).split('.')[0]
            # Path dictionary for the reference files
            _ref_path_dict_[file_name_gene] = ref_file_path
            ref_seq_count, _ref_length_dict_[file_name_gene] = get_reads_info(ref_file_path)
            # Calculate the average length
            _ref_length_dict_[file_name_gene] = int(_ref_length_dict_[file_name_gene] / max(ref_seq_count, 1))
    return _ref_length_dict_, _ref_path_dict_


# build a Hash table and store the reference in it
def make_ref_kmer_dict(reference_path: str, _kmer_size_: int, _ref_reverse_complement_: bool = False,
                       _pos_: bool = True, _print_: bool = True, ref_number: int = None, random_seed: int = 10) -> dict:
    files_list = get_file_list(reference_path)
    if not files_list:
        print('Reference is invalid. Please check the path!')
        sys.exit(1)
    ref_kmer_dict = defaultdict(list)
    gene_number_count = 0
    for file in files_list:
        file_name_gene = os.path.basename(file).split(".")[0]
        seq_total, need_to_read, ref_count = 0, None, 0
        if ref_number is not None:
            # Get the total number of sequences in the file
            with open(file, 'r', encoding='utf-8', errors='ignore') as f:
                seq_total = f.read().count('>')
            if seq_total <= ref_number:
                need_to_read = list(range(1, seq_total + 1))
            else:
                random.seed(random_seed)
                need_to_read = random.sample(list(range(1, seq_total + 1)), ref_number)
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        # Skip the first line >
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            gene_number_count += 1
            ref_count += 1
            if ref_number is not None:
                if ref_count not in need_to_read:
                    continue
                else:
                    need_to_read.remove(ref_count)
                if len(need_to_read) == 0:
                    break
            # Preserve the numbers and letters in the seq sequence Remove line breaks
            ref_seq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            for j in range(0, len(ref_seq) - _kmer_size_ + 1):
                kmer_info_list, ref_kmer = [], ref_seq[j:j + _kmer_size_]
                if ref_kmer in ref_kmer_dict:
                    kmer_info_list = ref_kmer_dict[ref_kmer]
                    # temp_list[0] is the number of occurrences of ref_kmer
                    kmer_info_list[0] += 1
                    if _pos_:
                        # temp_list[1] is the sum of the percentage positions where ref_kmer appears
                        kmer_info_list[1] += (j + 1) / len(ref_seq)
                    if file_name_gene not in kmer_info_list:
                        kmer_info_list.append(file_name_gene)
                else:
                    kmer_info_list = [1, (j + 1) / len(ref_seq), file_name_gene] if _pos_ else [1, 0,
                                                                                                file_name_gene]
                    # ref_kmer_dict: key->kmer_seq value->list
                    # list[0] Number of occurrences of kmer
                    # list[1] the sum of the percentage positions of kmer appearing on a reference gene
                    # list[2:] the genes in which the kmer appears
                    ref_kmer_dict[sys.intern(ref_kmer)] = kmer_info_list
            # _ref_reverse_complement_ reverse and complement the reference sequence
            if _ref_reverse_complement_:
                ref_seq = reverse_complement_limit(ref_seq)
                for j in range(0, len(ref_seq) - _kmer_size_ + 1):
                    kmer_info_list, ref_kmer = [], ref_seq[j:j + _kmer_size_]
                    if ref_kmer in ref_kmer_dict:
                        kmer_info_list = ref_kmer_dict[ref_kmer]
                        kmer_info_list[0] += 1
                        if _pos_:
                            # Negative value of pos for reverse sequence
                            kmer_info_list[1] -= (j + 1) / len(ref_seq)
                        if file_name_gene not in kmer_info_list:
                            kmer_info_list.append(file_name_gene)
                    else:
                        kmer_info_list = [1, -(j + 1) / len(ref_seq), file_name_gene] if _pos_ else [1, 0,
                                                                                                     file_name_gene]
                        ref_kmer_dict[sys.intern(ref_kmer)] = kmer_info_list
            if _print_:
                print('Mem.:', round(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, 2),
                      'G, Num. of Seq:', gene_number_count, end='\r')
            else:
                if gene_number_count % 1000 == 0:
                    print("INFO: Mem. used: {:.2f} G, Num. of Seq: {}".format(
                        psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, gene_number_count))
        infile.close()
    return ref_kmer_dict


def dynamic_limit_v2(_read_kmer_dict_: dict, smoother_level: int = 4, smoother_size: int = 64):
    count_list = [0] * smoother_size
    for x in _read_kmer_dict_:
        if _read_kmer_dict_[x][0] <= smoother_size:
            count_list[_read_kmer_dict_[x][0] - 1] += 1
    for i in range(smoother_level):
        for x in range(1, smoother_size - smoother_level + i):
            if count_list[x + smoother_level - 1 - i] < count_list[x - 1] \
                    and count_list[x] < count_list[x + smoother_level - i]:
                return x + 1
    return 2


# dynamic limit
def dynamic_limit(_read_kmer_dict_: dict, gen_avg_len: int, ref_rate: int = 1.5, list_size: int = 256) -> int:
    count_list = [0] * list_size
    for _read_kmer_ in _read_kmer_dict_:
        # _read_kmer_dict_[x][0] is the number of occurrences of kmer
        # the total number of all kmers with n occurrences is calculated
        if _read_kmer_dict_[_read_kmer_][0] <= list_size:
            count_list[_read_kmer_dict_[_read_kmer_][0] - 1] += 1
    F0, sum_f = len(_read_kmer_dict_), 0
    for x in range(list_size):
        sum_f += count_list[x]
        if (F0 - sum_f) / 2 < gen_avg_len * ref_rate:
            return max(2, x)
    return 2
