#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:59 下午
# @Author     : zzhen
# @File       : filter.py
# @Software   : PyCharm
# @Description: the file to implement the filter function
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import argparse
import glob
import gzip
import logging
import os
import platform
import sys
import time
import multiprocessing
from collections import defaultdict
from collections import deque
from itertools import chain
from sys import getsizeof, stderr
import shutil
import psutil

try:
    from reprlib import repr
except ImportError:
    pass

logger = logging.getLogger("easy353.filter")


def get_deep_sizeof(o, handlers=None, verbose=False):
    '''
    Returns the approximate memory footprint an object and all of its contents.
    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:
    handlers = {SomeContainerClass: iter, OtherContainerClass: OtherContainerClass.get_elements}
    '''
    if handlers is None:
        handlers = {}
    DEFAULT_HANDLERS = {tuple: iter, list: iter, deque: iter, dict: lambda d: chain.from_iterable(d.items()), set: iter,
                        frozenset: iter, }
    all_handlers = {**DEFAULT_HANDLERS, **handlers}  # user handlers take precedence
    seen = set()  # track which object id's have already been seen
    default_size = getsizeof(0)  # estimate sizeof object without __sizeof__

    def sizeof(obj):
        if id(obj) in seen:  # do not double count the same object
            return 0
        seen.add(id(obj))
        s = getsizeof(obj, default_size)
        if verbose:
            print(s, type(obj), repr(obj), file=stderr)
        for container_typ, container_handler in all_handlers.items():
            if isinstance(obj, container_typ):
                s += sum(map(sizeof, container_handler(obj)))
                break
        return s

    return sizeof(o)


# convert a DNA sequence to integer
def dnaseq2int(dna_seq: str, reverse_comp: bool = False):
    if reverse_comp:
        dna_seq = dna_seq.upper().translate(str.maketrans('ACGTU', '32100', 'RYMKSWHBVDN\n'))[::-1]
        return int(dna_seq, 4) if dna_seq else 0
    else:
        dna_seq = dna_seq.upper().translate(str.maketrans('ACGTU', '01233', 'RYMKSWHBVDN\n'))
        # seq is empty or composed of other characters
        if not dna_seq:
            return 0, 0
    return int(dna_seq, 4), len(dna_seq)


# convert bin seq to dna seq
def binseq2dna(bin_seq, seq_len):
    # 如果是二进制字符串，用int进行转换成整数
    # 如果是整数，直接当作10进制数继续处理
    decimal_num = int(bin_seq, 2) if isinstance(bin_seq, str) else bin_seq
    # 暗含了字符对应表 ('0123', 'ACGT')
    digits = 'ACGT'
    res = []
    while decimal_num:
        res.append(digits[decimal_num % 4])
        decimal_num >>= 2
    # 处理解码后的序列长度不够的情况
    return ''.join(res)[::-1].rjust(seq_len, 'A')


# Used to determine the file suffix name, used to select the specific read function
def biofile_type(path):
    suffix_dict = {'.gz': 0, '.fq': 1, '.fastq': 1, 'fa': 2, '.fas': 2, '.fasta': 2}
    # 设定默认值是3
    return suffix_dict.get(os.path.splitext(path)[-1].lower(), 3)


# if the path is a file,the file list is [path]
# if the path is a dir,the file list is [files under the path]
def get_file_lst(path: str) -> list:
    extension = ('.fasta', '.fas', '.fa', '.fna', '.ffn', '.frn', '.faa', '.fna')
    if os.path.isdir(path):
        return [_ for _ in glob.glob(os.path.join(path, '*')) if os.path.isfile(_) and _.endswith(extension)]
    elif os.path.isfile(path) and path.endswith(extension):
        return [path]
    else:
        logging.error('The path of ref is invalid!')
        sys.exit(-1)


def make_ref_kmer_dict(ref_kmer_dict: defaultdict, ref_path: str, k_size: int):
    # all the ref files
    # 对ref_path_lst进行sort 保证顺序
    ref_path_lst = sorted(get_file_lst(ref_path))
    seq_num_count = 0
    # 设置kmer_mask
    kmer_mask = (1 << (k_size << 1)) - 1
    # 获取参考文件列表长度
    ref_path_lst_len = len(ref_path_lst)
    for ref_index in range(0, ref_path_lst_len):
        ref_file = open(ref_path_lst[ref_index], 'rt')
        # 取出文件第一行
        next(ref_file)
        for line in ref_file:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = next(ref_file, None)
            # Preserve the numbers and letters in the seq sequence
            # Remove line breaks
            ref_seq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            # 基因数量+1
            seq_num_count += 1
            ref_bin_seq, seq_len = dnaseq2int(ref_seq)
            # 对ref序列进行反向互补
            ref_bin_seqs_lst = [ref_bin_seq, dnaseq2int(ref_seq, True)]
            for rc_flag, ref_bin_seq in enumerate(ref_bin_seqs_lst):
                bin_len = 22 + ref_path_lst_len
                for seq_index in range(seq_len - k_size + 1):
                    # 从前向后算
                    # 正反链信息 occurrence 文件信息
                    # 第1位：辅助位 第2位：正反链 第3-22位：occurrence 第23-(22+len(ref_path_lst)): 文件信息
                    bin_kmer_info = (1 << (bin_len - 1)) + (rc_flag << (bin_len - 2))
                    # 从ref_bin_seq中获取指定长度的kmer
                    bin_kmer = (ref_bin_seq >> (seq_index << 1)) & kmer_mask
                    if bin_kmer in ref_kmer_dict:
                        bin_kmer_info = ref_kmer_dict[bin_kmer]
                    # occurrence +1
                    bin_kmer_info += (1 << ref_path_lst_len)
                    # 设置文件位
                    bin_kmer_info |= (1 << (ref_path_lst_len - ref_index - 1))
                    ref_kmer_dict[bin_kmer] = bin_kmer_info
        logger.info('Num. of Hashed Ref Seq: {}'.format(seq_num_count))


# 根据read拆分成的kmer 判断该read可能属于的gene文件
# 返回一个文件index列表
def filter_read_for_ref(ref_kmer_dict: defaultdict, k_size: int, s_size: int, read_seqs_lst: list, ref_file_nums: int,
                        filter_pair_read: bool = False):
    if not read_seqs_lst:
        logger.warning('The read is empty!')
        return []
    # Initialize bit vector for storing matching file indices
    file_bits = 0
    paired_reads = len(read_seqs_lst) > 1
    # Convert read sequence to integer
    bin_read_seq, read_len = dnaseq2int(read_seqs_lst[0])
    # 当传入的序列是两条时，且需要处理filter_pair_read时
    bin_read_seqs_lst = [bin_read_seq, dnaseq2int(read_seqs_lst[1], True)] if filter_pair_read and paired_reads else [
        bin_read_seq]
    # kmer mask for extracting k-mers from read sequence
    kmer_mask = (1 << (k_size << 1)) - 1
    for bin_read_seq in bin_read_seqs_lst:
        for read_index in range(0, read_len - k_size - s_size + 2, s_size):
            read_kmer = (bin_read_seq >> (read_index << 1)) & kmer_mask
            if read_kmer in ref_kmer_dict:
                # 合并所有文件位置
                file_bits |= ref_kmer_dict[read_kmer]
    # Convert bit vector to list of file indices
    # 文件顺序是从左向右
    file_indices = [i for i in range(ref_file_nums) if (file_bits >> (ref_file_nums - 1 - i)) & 1]
    return file_indices


def reads_bin_filter(ref_path: str, ref_kmer_dict: defaultdict, k_size: int, s_size: int, fq_1: str, fq_2: str,
                     out_dir: str, p_id: int, p_count: int, filter_pair_read: bool = True):
    t_start, t_end, reads_count, ref_path_lst, file_dict = time.perf_counter(), 0, 0, sorted(get_file_lst(ref_path)), {}
    # detect whether the input is paired-end reads
    paired_reads = not fq_1 == fq_2
    # 设置输出文件
    for file_index, one_ref_path in enumerate(ref_path_lst):
        ref_file_name = os.path.splitext(os.path.basename(one_ref_path))[0]
        # out_fmt=1 输出是fastq out_fmt=0 输出是fasta
        # 将输出文件的句柄存放在dict中
        if not paired_reads:
            file_dict[file_index] = open(os.path.join(out_dir, ref_file_name + '.' + str(p_id + 1) + '.fasta'), 'w')
        else:
            # paired-end reads: the reads are written to two files
            file_dict[file_index] = open(os.path.join(out_dir, ref_file_name + '_R1' + '.' + str(p_id + 1) + '.fasta'),
                                         'w'), open(
                os.path.join(out_dir, ref_file_name + '_R2' + '.' + str(p_id + 1) + '.fasta'), 'w')
    # 判断文件是否是gz文件
    gz_file = True if fq_1.endswith('.gz') else False
    infile_1 = gzip.open(fq_1, 'rt') if gz_file else open(fq_1, 'rt')
    if paired_reads: infile_2 = gzip.open(fq_2, 'rt') if gz_file else open(fq_2, 'rt')
    for _ in infile_1:
        reads_count += 1
        tmp_rec1 = [_, next(infile_1, None), next(infile_1, None), next(infile_1, None)]
        if paired_reads:  tmp_rec2 = [next(infile_2, None), next(infile_2, None), next(infile_2, None),
                                      next(infile_2, None)]
        if reads_count % p_count == p_id:
            # write kmers into the files
            read_seqs_lst = [tmp_rec1[1], tmp_rec2[1]] if paired_reads and filter_pair_read else [tmp_rec1[1]]
            for file_index in filter_read_for_ref(ref_kmer_dict, k_size, s_size, read_seqs_lst, len(ref_path_lst),
                                                  filter_pair_read):
                if paired_reads:
                    file_dict[file_index][0].writelines(['>', tmp_rec1[0], tmp_rec1[1]])
                    file_dict[file_index][1].writelines(['>', tmp_rec2[0], tmp_rec2[1]])
                else:
                    file_dict[file_index].writelines(['>', tmp_rec1[0], tmp_rec1[1]])
            if (reads_count << int(paired_reads)) % 1000000 == 0:
                t_end = time.perf_counter()
                # 交换时间
                t_start, t_end = t_end, t_end - t_start
                logger.info(
                    'Handled {} M reads. {:.2f} s/M reads'.format((reads_count << int(paired_reads)) // 1000000,
                                                                  t_end))
    # 关闭文件句柄
    for _ in file_dict.keys():
        if paired_reads:
            file_dict[_][0].close()
            file_dict[_][1].close()
        else:
            file_dict[_].close()
    infile_1.close()
    if paired_reads: infile_2.close()


# logger = logging.getLogger(__file__)
# logger.setLevel(logging.DEBUG)
#
# # 建立一个filehandler来把日志记录在文件里，级别为debug以上
# fh = logging.FileHandler('filter.log')
# fh.setLevel(logging.DEBUG)
# fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(lineno)s : %(message)s',
#                                   datefmt='%Y-%m-%d %H:%M:%S'))
# # 将相应的handler添加在logger对象中
# logger.addHandler(fh)
# # 建立一个stream-handler来把日志打在CMD窗口上，级别为error以上
# ch = logging.StreamHandler(sys.stdout)
# # 当使用silent mode时，只输出error日志
# print_level = logging.DEBUG
# ch.setLevel(print_level)
# ch.setFormatter(logging.Formatter('%(message)s'))
# logger.addHandler(ch)
#
# ref_path = '/Users/zzhen/Desktop/test/Glycine_353'
# fq1 = '/Users/zzhen/Desktop/test/Gmax_sim_1.fastq.gz'
# fq2 = '/Users/zzhen/Desktop/test/Gmax_sim_2.fastq.gz'
# out_dir = '/Users/zzhen/Desktop/test/my_test'
# ref_kmer_dict = defaultdict(int)
#
# t_start = time.perf_counter()
# make_ref_kmer_dict(ref_kmer_dict, ref_path, 21, True, False)
# reads_bin_filter(ref_path, ref_kmer_dict, 21, 1, fq1, fq2, out_dir, 0, 1, False)
#
# t_end = time.perf_counter()
#
# logger.info('Time used for filter: {:.2f} s'.format(t_end - t_start))

if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage='%(prog)s [options]',
                                   description="Easy353's filter zzhen@sculab")
    # pars.add_argument('-fq1', dest='fq_file_1', type=str, help='Input sequencing file(s) -1 (*.fq/.gz).',
    #                   required=False)
    # pars.add_argument('-fq2', dest='fq_file_2', type=str, help='Input sequencing file(s) -2 (*.fq/.gz).',
    #                   required=False)
    # pars.add_argument('-r', dest='reference', type=str, help='Input a ref file or directory.')
    # pars.add_argument('-o', dest='out_dir', type=str, help='Output directory.', default='filter_out')
    # pars.add_argument('-fk', dest='filter_kmer', type=int, help='Kmer setting for filtering reads. Default:31',
    #                   default=31)
    # pars.add_argument('-s', dest='step_length', type=int,
    #                   help='Length of the sliding window on the reads. Default:1', default=1)
    # pars.add_argument('-ft', dest='filter_thread', type=int,
    #                   help='Threads setting for filtering reads. Default:1', default=1)
    # pars.add_argument('-fpr', dest='filter_pair_read', action='store_true',
    #                   help='fpr mode:get more filtered reads. If yes (more reads) or not (less time).')
    # pars.add_argument('-log', dest='log_file', type=str, help='The log file.', default='filter.log')
    # pars.add_argument('-silent', dest='silent', action='store_true',
    #                   help='If using silent mode, there will not output any info!')
    args = pars.parse_args()

    args.fq_file_1 = '/Volumes/zzhen/data/work_data/Easy353/empirical_data/tran_apiaceae/Anthriscus_sylvestris.1.fq.gz'
    args.fq_file_2 = '/Volumes/zzhen/data/work_data/Easy353/empirical_data/tran_apiaceae/Anthriscus_sylvestris.2.fq.gz'
    args.reference = '/Users/zzhen/Desktop/test/Apiaceae353'
    args.out_dir = '/Users/zzhen/Desktop/test/filter_out'
    args.filter_kmer = 21
    args.step_length = 1
    args.filter_thread = 8
    args.filter_pair_read = False
    args.log_file = 'filter.log'
    args.silent = False
    # print(args)

    # print(get_file_lst(args.reference))
    # print(sorted(get_file_lst(args.reference)))
    # for i, j in enumerate(sorted(get_file_lst(args.reference))):
    #     print(i, j)

    # set the logger
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.DEBUG)

    # 建立一个filehandler来把日志记录在文件里，级别为debug以上
    fh = logging.FileHandler(args.log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(lineno)s : %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S'))
    # 将相应的handler添加在logger对象中
    logger.addHandler(fh)
    # 建立一个stream-handler来把日志打在CMD窗口上，级别为error以上
    ch = logging.StreamHandler(sys.stdout)
    # 当使用silent mode时，只输出error日志
    print_level = logging.ERROR if args.silent else logging.DEBUG
    ch.setLevel(print_level)
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)

    #  set the start method of multiprocessing
    if platform.system() in ('Linux', 'Darwin'):
        # fork方式 子进程可以共享父进程的资源
        multiprocessing.set_start_method('fork')
    else:
        # spawn方式 子进程会copy一份父进程的资源，导致内存占用较大
        multiprocessing.set_start_method('spawn')

    # 先检查文件夹是否都存在
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    else:
        #    输出文件夹不允许存在
        logger.warning('The output directory: {} exists, it will be deleted and remade!'.format(args.out_dir))
        shutil.rmtree(args.out_dir)
        os.makedirs(args.out_dir)

    logger.info('Getting information from references')

    t_hash_start = time.perf_counter()
    ref_kmer_dict = defaultdict(int)
    # build hash table
    logger.info('Building binary hash table')

    make_ref_kmer_dict(ref_kmer_dict, args.reference, args.filter_kmer)
    logger.info('Hash dictionary has been made')
    t_hash_end = time.perf_counter()
    logger.info('Time used for building hash table: {:.2f} s'.format(t_hash_end - t_hash_start))

    logger.info('The memory usage of the current process is: {:.2f} GB'.
                format(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024))
    logger.info('Filter reads from fq_files based on hash table with {} thread'.format(args.filter_thread))

    # 处理输入的测序文件
    if not args.fq_file_1 and not args.fq_file_2:
        logger.error('There were no sequencing files inputted.')
        sys.exit(-1)
    elif args.fq_file_1 and not args.fq_file_2:
        args.fq_file_2 = args.fq_file_1
    elif args.fq_file_2 and not args.fq_file_1:
        args.fq_file_1 = args.fq_file_2
    if args.filter_thread == 1:
        # when filter_thread is 1
        reads_bin_filter(args.reference, ref_kmer_dict, args.filter_kmer, args.step_length,
                         args.fq_file_1, args.fq_file_2,
                         args.out_dir, 0, args.filter_thread, args.filter_pair_read)
    else:
        process_list = []
        for p_id in range(args.filter_thread):
            p = multiprocessing.Process(target=reads_bin_filter,
                                        args=(args.reference, ref_kmer_dict, args.filter_kmer, args.step_length,
                                              args.fq_file_1, args.fq_file_2,
                                              args.out_dir, p_id, args.filter_thread, args.filter_pair_read,))
            process_list.append(p)
        for p in process_list: p.start()
        for p in process_list: p.join()

    # 合并多进程产生的文件
    ref_path_lst = get_file_lst(args.reference)
    paired_reads = True if args.fq_file_1 != args.fq_file_2 else False
    for ref_file in ref_path_lst:
        gene_name = os.path.splitext(os.path.basename(ref_file))[0]
        prefixes = [gene_name + "_R1", gene_name + "_R2"] if paired_reads else [gene_name]
        for prefix in prefixes:
            with open(os.path.join(args.out_dir, prefix + ".fasta"), 'wt') as merge_file:
                for file_name in [os.path.join(args.out_dir, prefix + "." + str(i + 1) + ".fasta") for i in
                                  range(args.filter_thread)]:
                    with open(file_name, 'rt') as infile:
                        shutil.copyfileobj(infile, merge_file)
                    os.remove(file_name)

    t_filter_end = time.perf_counter()
    logger.info('Time used for filter: {:.2f} s'.format(t_filter_end - t_hash_end))
    logger.info('Total time used for filter: {:.2f} s'.format(t_filter_end - t_hash_start))
