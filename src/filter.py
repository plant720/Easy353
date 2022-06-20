#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:59 下午
# @Author     : zzhen
# @File       : filter_1.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import collections
import gzip
import multiprocessing
import os
import shutil
import sys
import time
import psutil
import gc
import argparse

# import my python file
from src.utils import bytes_str, reverse_complement_limit, make_ref_kmer_dict, get_file_size, get_ref_info, total_size, \
    get_reads_info,log


# 根据read_seq判断该read产生的kmer包含在哪些参考文件中,返回文件名/基因名 即os.path.basename(file_path).split('.')[0]
# ref_hash_dict是参考序列的dict
# 如果建立_ref_hash_dict_时,对参考序列进行了反向互补，则在过滤read时不需要;反之同理
def filter_read_for_gene(_ref_kmer_dict_: dict, _kmer_size_: int, _step_size_: int, _read_seq_: str,
                         _read_reverse_complement_: bool = True) -> set:
    ref_gene_name_list = []
    # 由于换行符的存在,不需要len-kmer_size+1
    for j in range(0, len(_read_seq_) - _kmer_size_, _step_size_):
        read_kmer = _read_seq_[j:j + _kmer_size_]
        if read_kmer in _ref_kmer_dict_:
            # 将该kmer存在的参考文件添加到tmplist
            ref_gene_name_list.extend(_ref_kmer_dict_.get(read_kmer)[2:])
    if _read_reverse_complement_:
        _read_seq_ = reverse_complement_limit(_read_seq_)
        for j in range(0, len(_read_seq_) - _kmer_size_, _step_size_):
            read_kmer = _read_seq_[j:j + _kmer_size_]
            if read_kmer in _ref_kmer_dict_:
                ref_gene_name_list.extend(_ref_kmer_dict_.get(read_kmer)[2:])
    return set(ref_gene_name_list)


# _ref_hash_dict_是参考序列构建的哈希字典
# _step_size_ the length of the sliding window on the reads
# 如果建立_ref_hash_dict_时,对参考序列进行了反向互补，则在过滤read时不需要;反之同理
def do_reads_filter(_ref_kmer_dict_: dict, _kmer_size_: int, _step_size_: int, _file_path_: str, _out_dir_: str,
                    _read_reverse_complement_, _print_=True) -> None:
    _time_start_, _time_end_, reads_count = time.time(), 0, 0
    # 判断文件是否是gz文件
    bytes_type = _file_path_[-3:].lower() == ".gz"
    infile_stream = gzip.open(_file_path_, 'rb') if bytes_type else open(_file_path_, 'rb')
    for _ in infile_stream:
        reads_count += 1
        temp_rec = [_, infile_stream.readline(), infile_stream.readline(), infile_stream.readline()]
        for file_name in filter_read_for_gene(_ref_kmer_dict_, _kmer_size_, _step_size_,
                                              bytes_str(temp_rec[1]),
                                              _read_reverse_complement_):
            with open(os.path.join(_out_dir_, file_name + ".fasta"), "a+") as outfile:
                outfile.writelines(['>', bytes_str(temp_rec[0]), bytes_str(temp_rec[1])])
        if reads_count % 1000000 == 0:
            _time_end_ = time.time()
            _time_start_, _time_end_ = _time_end_, _time_end_ - _time_start_
            if _print_:
                print('handled\t', reads_count // 1000000, 'm reads, ', round(_time_end_, 2), 's/m reads', sep="",
                      end='\r')
    infile_stream.close()


# 如果建立_ref_hash_dict_时,对参考序列进行了反向互补，则在过滤read时不需要;反之同理
# 出于内存考虑默认设置ref_reverse_complement=False 故_read_reverse_complement_=True
# t_id and t_count 适用于处理多进程读文件
def do_pair_reads_filter(_ref_kmer_dict_: dict, _kmer_size_: int, _step_size_: int, fq_file_1: str,
                         fq_file_2: str, _out_dir_: str, t_id: int, t_count: int,
                         _read_reverse_complement_, _print_=True) -> None:
    _time_start_, _time_end_, reads_count = time.time(), 0, 0
    bytes_type = fq_file_1[-3:].lower() == ".gz"
    if fq_file_1 == fq_file_2:
        print("Error: the paired reads file is same!")
        exit(1)
    infile_1 = gzip.open(fq_file_1, 'rb') if bytes_type else open(fq_file_1, 'rb')
    infile_2 = gzip.open(fq_file_2, 'rb') if bytes_type else open(fq_file_2, 'rb')
    for _ in infile_1:
        reads_count += 1
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
        if reads_count % t_count == t_id:
            # 对于该kmer存在的所有参考文件 进行文件写入
            for file_name in filter_read_for_gene(_ref_kmer_dict_, _kmer_size_, _step_size_,
                                                  bytes_str(temp_rec1[1]),
                                                  _read_reverse_complement_):
                with open(os.path.join(_out_dir_, file_name + ".fasta"), "a+") as _out_file:
                    _out_file.writelines(['>', bytes_str(temp_rec1[0]), bytes_str(temp_rec1[1])])
                    _out_file.writelines(['>', bytes_str(temp_rec2[0]), bytes_str(temp_rec2[1])])
        if reads_count * 2 % 1000000 == 0:
            _time_end_ = time.time()
            _time_start_, _time_end_ = _time_end_, _time_end_ - _time_start_
            if _print_:
                print('handled\t', reads_count * 2 // 1000000, 'm reads, ', round(_time_end_, 2), 's/m reads',
                      sep="", end='\r')
    infile_1.close()
    infile_2.close()


# 用于对过滤后的数据较大的情况增大kmer进行过滤
def re_filter_reads(_ref_kmer_dict_: dict, _kmer_size_: int, _step_size_: int, gene_name: str, _out_dir_: str,
                    _read_reverse_complement_=True) -> tuple:
    refilter_reads_count = 0
    refilter_reads_length = 0
    out_file_stream = open(os.path.join(_out_dir_, gene_name + ".fasta"), "w+")
    infile_stream = open(os.path.join(_out_dir_, "big_reads", gene_name + ".fasta"), 'r')
    for _ in infile_stream:
        read_seq = infile_stream.readline()
        filtered_flag = False
        for j in range(0, len(read_seq) - _kmer_size_, _step_size_):
            kmer = read_seq[j:j + _kmer_size_]
            if kmer in _ref_kmer_dict_:
                filtered_flag = True
                break
        # 当该条序列没有过滤成功并且需要处理read反向互补的情况时
        if not filtered_flag and _read_reverse_complement_:
            temp_seq = reverse_complement_limit(read_seq)
            # 因为有一个换行符 所以不需要-1
            for j in range(0, len(temp_seq) - _kmer_size_, _step_size_):
                kmer = temp_seq[j:j + _kmer_size_]
                if kmer in _ref_kmer_dict_:
                    filtered_flag = True
                    break
        if filtered_flag:
            refilter_reads_count += 1
            refilter_reads_length += len(read_seq)
            out_file_stream.writelines([_, read_seq])
    out_file_stream.close()
    infile_stream.close()
    return refilter_reads_count, refilter_reads_length


# 对单个基因文件进行过滤
def refilter_one_gene(gene_name: str, gene_avg_len: int, _refilter_kmer_size_: int,
                      _reference_path_: str, _out_dir_: str, _refilter_step_size_: int = 1):
    _refilter_gene_info_dict_ = collections.defaultdict(dict)
    _refilter_gene_info_dict_[gene_name]["ref_avg_length"] = gene_avg_len
    while True:
        print("re-filter:", gene_name, 'with k =', _refilter_kmer_size_)
        # 对参考序列进行reverse_complement 则不对read进行反向互补
        ref_kmer_dict = make_ref_kmer_dict(_reference_path_, _refilter_kmer_size_, _ref_reverse_complement_=True,
                                           _print_=False)
        filter_reads_count, filter_reads_length = re_filter_reads(ref_kmer_dict, _refilter_kmer_size_,
                                                                  _refilter_step_size_, gene_name,
                                                                  _out_dir_, _read_reverse_complement_=False)

        _refilter_gene_info_dict_[gene_name]["filter_reads_count"] = filter_reads_count
        _refilter_gene_info_dict_[gene_name]["filter_reads_length"] = filter_reads_length
        _refilter_gene_info_dict_[gene_name]["filter_kmer"] = _refilter_kmer_size_

        if filter_reads_length / gene_avg_len < 512 or get_file_size(
                os.path.join(_out_dir_, gene_name + ".fasta")) < 8 or _refilter_kmer_size_ >= 75:
            break
        else:
            # 将之前big_reads下的文件删除，然后将过滤后的文件移动到big_reads下
            if os.path.isfile(os.path.join(_out_dir_, "big_reads", gene_name + ".fasta")):
                os.remove(os.path.join(_out_dir_, "big_reads", gene_name + ".fasta"))
            shutil.move(os.path.join(_out_dir_, gene_name + ".fasta"),
                        os.path.join(_out_dir_, 'big_reads', gene_name + ".fasta"))
            _refilter_kmer_size_ += 2
    return _refilter_gene_info_dict_


# _clear_参数用于清理输出文件夹下的文件
def filter_flow(_read_data_tuple_: tuple, _out_dir_: str, _reference_path_: str, _kmer_size_: int,
                _step_size_: int = 1, _ref_reverse_complement_: bool = False,
                _read_reverse_complement_: bool = False, refilter: bool = True,
                _clear_: bool = True, _pos_: bool = True,
                _paired_reads_: bool = True, _thread_for_filter_: int = 4, _print_: bool = True):
    # 初始化操作
    if not os.path.isdir(_out_dir_):
        os.makedirs(_out_dir_)
    print("Getting information from references")
    # 当参考序列不进行反向互补时，read进反向互补
    if not _ref_reverse_complement_:
        _read_reverse_complement_ = True

    _time_build_hash_start_ = time.time()
    # ref_length_dict 记录reference文件的read平均长度
    # ref_path_dict 记录reference文件的路径信息
    ref_length_dict, ref_path_dict = get_ref_info(_reference_path_)
    # 构建哈西字典
    print("Building hash table")

    # 处理reference 并生成具体kmer存储在ref_kmer_dict
    # 用于记录记录所有的reference信息  字典的key是kmer seq value是一个list,
    # list[0] 是kmer出现次数
    # list[1] 是kmer出现在某条参考基因上的百分比位置之和
    # list[2:]是该kmer在哪些文件(基因)中出现
    ref_kmer_dict = make_ref_kmer_dict(_reference_path_, _kmer_size_, _ref_reverse_complement_, _pos_, _print_)
    print("Hash dictionary has been made")
    _time_build_hash_end = time.time()
    print('Time used for building hash table: {} s'.format(round(_time_build_hash_end - _time_build_hash_start_, 2)))

    # _clear_用于将输出文件夹中的原有同名文件删除,新建文件
    if _clear_:
        for key in ref_path_dict:
            with open(os.path.join(_out_dir_, key + ".fasta"), 'w'):
                pass
    print('Filter reads from fq_files based on hash table')

    _time_filter_start_ = time.time()
    # 当测序数据是unpaired_end时
    if not _paired_reads_:
        _unpaired_file_list_ = _read_data_tuple_[0]
        # 当建立的hash表的内存超过系统内存一半时，使用单线程进行过滤
        if _thread_for_filter_ == 1 or total_size(ref_kmer_dict) / 1024 / 1024 / 1024 * _thread_for_filter_ > \
                psutil.virtual_memory().total / 1024 / 1024 / 1024 * 0.5:
            for unpaired_file in _unpaired_file_list_:
                do_reads_filter(ref_kmer_dict, _kmer_size_, _step_size_, unpaired_file,
                                _out_dir_, _read_reverse_complement_, _print_)
        else:
            # 创建进程池
            pool = multiprocessing.Pool(min(len(_unpaired_file_list_), _thread_for_filter_))
            # 创建进程池中的进程
            for unpaired_file in _unpaired_file_list_:
                pool.apply_async(do_reads_filter, args=(ref_kmer_dict, _kmer_size_, _step_size_, unpaired_file,
                                                        _out_dir_, _read_reverse_complement_, _print_), )
            pool.close()
            pool.join()
    # 当测序数据是paired-end时
    if _paired_reads_:
        paired_file_list_one = _read_data_tuple_[0]
        paired_file_list_two = _read_data_tuple_[1]
        if len(paired_file_list_one) == len(paired_file_list_two):
            for i in range(len(paired_file_list_one)):
                # 处理根据内存占比，选用多进程或是单进程
                if _thread_for_filter_ == 1 or total_size(ref_kmer_dict) / 1024 / 1024 / 1024 * _thread_for_filter_ > \
                        psutil.virtual_memory().total / 1024 / 1024 / 1024 * 0.5:
                    # 当设定为单线程或
                    do_pair_reads_filter(ref_kmer_dict, _kmer_size_, _step_size_,
                                         paired_file_list_one[i], paired_file_list_two[i],
                                         _out_dir_, 0, 1, _read_reverse_complement_, _print_)
                else:
                    # 创建进程池
                    pool = multiprocessing.Pool(_thread_for_filter_)
                    # 创建进程池中的进程
                    # 多进程读取一个文件
                    for j in range(_thread_for_filter_):
                        pool.apply_async(func=do_pair_reads_filter,
                                         args=(ref_kmer_dict, _kmer_size_, _step_size_,
                                               paired_file_list_one[i], paired_file_list_two[i],
                                               _out_dir_, j, _thread_for_filter_, _read_reverse_complement_, _print_,))
                    pool.close()
                    pool.join()
        else:
            print("The number of paired_end fastq files is not equal")
            exit(0)
    _time_filter_end_ = time.time()
    print('Time used for filter: {} s'.format(round(_time_filter_end_ - _time_filter_start_, 2)))
    # 过滤完后 手动清理内存
    del ref_kmer_dict
    gc.collect()

    # 根据过滤后的数据 获取过滤出的reads数量和reads总长度
    filter_gene_info_dict = collections.defaultdict(dict)
    for gene_name, ref_avg_length in ref_length_dict.items():
        filter_gene_info_dict[gene_name]["ref_avg_length"] = ref_avg_length
        filter_gene_info_dict[gene_name]["filter_kmer"] = _kmer_size_
        filter_gene_info_dict[gene_name]["filter_reads_count"], filter_gene_info_dict[gene_name][
            "filter_reads_length"] = get_reads_info(os.path.join(_out_dir_, gene_name + ".fasta"))

    # 当需要re-filter时
    if refilter:
        # 记录需要过滤的基因 对应参考基因平均长度和参考序列所在文件 key:gene_name value:[avg_len,gene_path]
        refilter_gene_info_dict = collections.defaultdict(dict)
        if not os.path.isdir(os.path.join(_out_dir_, 'big_reads')):
            os.makedirs(os.path.join(_out_dir_, 'big_reads'))
        # ref_length_dict 记录reference文件的read平均长度
        for gene_name, gene_info in filter_gene_info_dict.items():
            if os.path.isfile(os.path.join(_out_dir_, gene_name + ".fasta")):
                # 如果覆盖深度>512x 或过滤后的文件>8m
                # 该处覆盖深度采用简单的 c = LN / G
                if (gene_info["filter_reads_length"] / gene_info.get("ref_avg_length") > 512) \
                        and get_file_size(os.path.join(_out_dir_, gene_name + ".fasta")) > 8:
                    # 将需要过滤的文件移动到一个文件夹中
                    shutil.move(os.path.join(_out_dir_, gene_name + ".fasta"),
                                os.path.join(_out_dir_, 'big_reads', gene_name + ".fasta"))
                    refilter_gene_info_dict[gene_name] = gene_info
        # 再次清理内存
        gc.collect()
        if not refilter_gene_info_dict:
            print("Big reads files don't exist, no need to reFilter!")
        else:
            print("reFilter")
            #  当需要进行重新过滤时,一般单个文件较小,可以考虑使用多进程处理重过滤
            # py3.8 对于 macOS,默认启动方式是spawn,m1上会出错
            # if sys.platform in ['darwin']:
            #     multiprocessing.set_start_method('fork')
            # 创建进程池
            pool = multiprocessing.Pool(processes=min(_thread_for_filter_, len(refilter_gene_info_dict)))
            # 创建进程池中的进程
            for gene_name, gene_info in refilter_gene_info_dict.items():
                pool.apply_async(func=refilter_one_gene,
                                 args=(gene_name, gene_info["ref_avg_length"],
                                       _kmer_size_ + 2, ref_path_dict[gene_name],
                                       _out_dir_, 1,),
                                 callback=refilter_gene_info_dict.update)
            pool.close()
            pool.join()
        # 输出log信息
        filter_gene_info_dict.update(refilter_gene_info_dict)
        log_file = os.path.join(_out_dir_, "filter_log.csv")
        log(log_file, "gene_id", "ref_avg_length", "filter_reads_count", "filter_kmer")
        for gene_name, gene_info in filter_gene_info_dict.items():
            log(log_file, gene_name, gene_info["ref_avg_length"], gene_info["filter_reads_count"],
                      gene_info["filter_kmer"])
        # 返回过滤后的基因信息
        return filter_gene_info_dict


if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage="%(prog)s [options]",
                                   description="Easy353 zzhen@sculab")
    pars.add_argument("-1", dest="fq_file_1", type=str, nargs="+",
                      help="Input file(s) with forward paired-end reads (*.fq/.gz/.tar.gz).", required=False)
    pars.add_argument("-2", dest="fq_file_2", type=str, nargs="+",
                      help="Input file(s) with reverse paired-end reads (*.fq/.gz/.tar.gz).", required=False)
    pars.add_argument("-u", dest="unpaired_fq_file", type=str,
                      help="Input file(s) with unpaired (single-end) reads.", required=False, nargs="+")
    pars.add_argument("-r", dest="reference", type=str, help="Input a file(directory) with references", required=True)
    pars.add_argument("-o", dest="output_dir", type=str, help="Output directory.", required=False,
                      default="easy353_output")
    pars.add_argument("-k", dest="filter_kmer", type=int, help="Kmer setting for filtering reads. Default:31",
                      default=31)
    pars.add_argument("-s", dest="step_length", type=int,
                      help="Step length of the sliding window on the reads. Default:1", default=1)
    pars.add_argument("-t", dest="filter_thread", type=int,
                      help="Threads setting for filtering reads. Default:4", default=4)
    pars.add_argument("-fast", dest="fast", action="store_true", help="Whether to use fast mode.")
    args = pars.parse_args()

    # python filter.py -1 /Users/zzhen/Desktop/ara.1.fq -2 /Users/zzhen/Desktop/ara.2.fq
    # -r /Users/zzhen/Desktop/ref -o /Users/zzhen/Desktop/result -k 31 -s 1 -t 4 -fast
    fastq_files = tuple()
    _paired_reads_ = False
    if args.unpaired_fq_file:
        fastq_files = (args.unpaired_fq_file,)
    if args.fq_file_1 and args.fq_file_2:
        fastq_files = (args.fq_file_1, args.fq_file_2)
        _paired_reads_ = True

    filter_flow(_read_data_tuple_=fastq_files, _out_dir_=args.output_dir,
                _reference_path_=args.reference, _kmer_size_=args.filter_kmer, _step_size_=args.step_length,
                _ref_reverse_complement_=args.fast, _paired_reads_=_paired_reads_,
                _thread_for_filter_=args.filter_thread)
