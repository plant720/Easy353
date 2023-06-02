# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/4/1 22:25
# @Author     : zzhen
# @File       : assemble.py
# @Software   : PyCharm
# @Description: the file to implement the assembly function
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.

import os
import sys
import time
import shutil
import logging
import platform
import argparse
import multiprocessing
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
# import the function implemented in the filter_yy.py file
from Easy353Lib.filter import dnaseq2int, binseq2dna, get_file_lst
# import functions implemented in the node.py file
from Easy353Lib.node import DeBruijnGraph

logger = logging.getLogger("easy353.assemble")


# import ctypes
#
# utils = ctypes.CDLL('./clib/utils_mac.so') if __name__ == '__main__' else ctypes.CDLL('Easy353Lib/clib/utils_mac.so')
# utils.optimized_lcs_len.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
# utils.optimized_lcs_len.restype = ctypes.c_int
#
#
# def get_ref_sequence(ref_file_path: str) -> str:
#     ref_sequence = ""
#     with open(ref_file_path, 'rt') as ref_file:
#         for line in ref_file:
#             if line.startswith(">"):
#                 continue
#             ref_sequence += line
#     return ref_sequence.upper().translate(str.maketrans("ACGTU\n", "ACGTT$", "RYMKSWHBVDN"))
#
#
# # 根据每一条read和ref序列的overlap长度 产生一个用于拼接的kmer值
# def get_dynamic_kmer(fa_1: str, fa_2: str, ref_file_path: str, min_overlap_len: int = 21,
#                      max_overlap_len: int = 51, max_read_num: int = 2000):
#     overlap_len_lst = []
#     # get the reference sequence
#     ref_sequence = get_ref_sequence(ref_file_path).encode("utf-8")
#     paired_read = not fa_1 == fa_2
#     # get the read sequence
#     read_file_1 = open(fa_1, 'rt')
#     if paired_read:
#         read_file_2 = open(fa_2, 'rt')
#     for line in read_file_1:
#         if line[0] == ">":
#             continue
#         if paired_read:
#             next(read_file_2)
#         reads_seq = line + "$" + next(read_file_2) if paired_read else line
#         reads_seq = reads_seq.encode("utf-8")
#         # get the overlap length
#         overlap_len_lst.append(utils.optimized_lcs_len(ref_sequence, reads_seq, min_overlap_len, max_overlap_len))
#         max_read_num -= 1
#         if max_read_num <= 0:
#             break
#     if sum(overlap_len_lst) / len(overlap_len_lst) < min_overlap_len:
#         return min_overlap_len
#     overlap_len_lst = [i for i in overlap_len_lst if i > 0]
#     return int(sum(overlap_len_lst) / len(overlap_len_lst))
#
#
# # # 从read_File中选定一定数量的reads，同参考序列进行overlap，用于判断参考序列同contig的方向
# # def get_ref_direction(fa_1: str, fa_2: str, ref_file_path: str, read_num: int = 200):
# #     # get the reference sequence
# #     ref_sequence = get_ref_sequence(ref_file_path).encode("utf-8")
# #     forward_overlap_len_sum, reverse_overlap_len_sum = 0, 0
# #     paired_read = not fa_1 == fa_2
# #     # get the read sequence
# #     forward_file = open(fa_1, 'rt')
# #     if paired_read: reverse_file = open(fa_2, 'rt')
# #     for forward_read in forward_file:
# #         if forward_read[0] == ">":
# #             continue
# #         if paired_read:
# #             next(reverse_file)
# #         reverse_read = next(reverse_file) if paired_read else dna_reverse_complement(forward_read)
# #         # get the overlap length
# #         forward_overlap_len_sum += utils.optimized_lcs_len(ref_sequence, forward_read.encode("utf-8"), 0, 500)
# #         reverse_overlap_len_sum += utils.optimized_lcs_len(ref_sequence, reverse_read.encode("utf-8"), 0, 500)
# #         read_num -= 1
# #         if read_num <= 0:
# #             break
# #     return False if forward_overlap_len_sum > reverse_overlap_len_sum else True


# # get the avg length of the reference seqs
# def get_avg_ref_len(ref_path: str) -> int:
#     ref_len_lst = []
#     with open(ref_path, 'rt') as ref_file:
#         for line in ref_file:
#             tmp_seq = []
#             while line and line[0] != '>':
#                 tmp_seq.append(line)
#                 line = next(ref_file, None)
#             # Only keep the letters in the ref sequence
#             ref_seq = ''.join(filter(str.isalpha, ''.join(tmp_seq).upper()))
#             ref_len_lst.append(len(ref_seq))
#     return int(sum(ref_len_lst) / len(ref_len_lst))

# return the reverse complementary sequence of a DNA sequence.
def dna_reverse_complement(_dna_: str) -> str:
    # use str.maketrans to build a translation table
    return _dna_.upper().translate(str.maketrans("ACGTU", "TGCAA", "RYMKSWHBVDN\n"))[::-1]


# extract all sequences from a fasta file and return a list of tuples (id,seq)
def extract_seq_from_fasta_file(ref_file_path: str) -> list:
    sequences = []
    with open(ref_file_path, 'rt') as ref_file:
        seq_id, seq = "", ""
        for line in ref_file:
            if line.startswith(">"):
                if seq != "":
                    sequences.append((seq_id, seq))
                    seq = ""
                seq_id = line[1:].strip()
                continue
            seq += line.strip()
        if seq != "":
            sequences.append((seq_id, seq))
    return sequences


# build a hash table based on the reference file which is a fasta file
def build_ref_hash_table(ref_kmer_dict: defaultdict, ref_file_path: str, k_size: int, get_reverse_comp: bool = False):
    # set the kmer mask
    kmer_mask = (1 << (k_size << 1)) - 1
    ref_file = open(ref_file_path, 'rt')
    # skip the first line
    next(ref_file)
    for line in ref_file:
        tmp_seq = []
        while line and line[0] != '>':
            tmp_seq.append(line)
            line = next(ref_file, None)
        # Only keep the letters in the ref sequence
        ref_seq = ''.join(filter(str.isalpha, ''.join(tmp_seq).upper()))
        # convert the DNA sequence to a binary sequence
        ref_bin_seq, seq_len = dnaseq2int(ref_seq)
        # add the reverse complementary sequence to the list
        ref_bin_seqs_lst = [ref_bin_seq, dnaseq2int(ref_seq, True)] if get_reverse_comp else [ref_bin_seq]
        for rc_flag, ref_bin_seq in enumerate(ref_bin_seqs_lst):
            for seq_index in range(seq_len - k_size + 1):
                # 从ref_bin_seq中获取指定长度的kmer
                bin_kmer = (ref_bin_seq >> (seq_index << 1)) & kmer_mask
                if bin_kmer in ref_kmer_dict:
                    ref_kmer_dict[bin_kmer] += 1
                else:
                    ref_kmer_dict[bin_kmer] = 1
    ref_file.close()
    # log information
    logger.debug(
        "The hash table based on the {} has been built.".format(os.path.splitext(os.path.basename(ref_file_path))[0]))


# build a hash table based on the read file which is a fasta file
# read_kmer_dict: key:read_kmer value:[read_count,read_pos]
def build_read_hash_table(fa_1: str, fa_2: str, read_kmer_dict: defaultdict, k_size: int):
    logger.debug("Building the hash table based on the read file...")
    t_start, t_end, kmer_mask = time.perf_counter(), 0, (1 << (k_size << 1)) - 1
    # detect whether the input is paired-end sequencing data
    paired_reads = not fa_1 == fa_2
    infile_1 = open(fa_1, "rt")
    if paired_reads: infile_2 = open(fa_2, "rt")
    # generate read sequences from the read file
    for _ in infile_1:
        # read_seq_info_lst is a list which contains the binary sequence and the length of the sequence
        infile_1_seq = next(infile_1)
        read_seqs_info_lst = [dnaseq2int(infile_1_seq), dnaseq2int(dna_reverse_complement(infile_1_seq))]
        if paired_reads:
            next(infile_2)
            # reverse complement the second read firstly and then convert it to a binary sequence and add it to the list
            # the paired reads are not completely overlapped
            # so the length of reverse_comp seq of second read is not equal to the length of the first read
            infile_2_seq = next(infile_2)
            read_seqs_info_lst.append(dnaseq2int(infile_2_seq))
            read_seqs_info_lst.append(dnaseq2int(dna_reverse_complement(infile_2_seq)))
        for read_seq, seq_len in read_seqs_info_lst:
            for seq_index in range(0, seq_len - k_size + 1):
                # get the kmer
                bin_kmer = (read_seq >> (2 * seq_index)) & kmer_mask
                # if bin_kmer in read_kmer_dict, the occurrence of the kmer is added 1, otherwise, the occurrence is set to 1
                if bin_kmer in read_kmer_dict:
                    read_kmer_dict[bin_kmer] += 1
                else:
                    read_kmer_dict[bin_kmer] = 1
    # close the file object
    infile_1.close()
    if paired_reads: infile_2.close()
    t_end = time.perf_counter()
    logger.debug("The read kmer dict has been built, used {:.2f}s".format(t_end - t_start))


def assemble_seq_from_dbg_v2(dbg_graph: DeBruijnGraph, ref_file: str, change_seed: int = 32,
                             max_iterations: int = 1000):
    ref_kmer_dict = defaultdict(int)
    k_size = dbg_graph.k_size
    best_seq, best_seed, seed_num = '', '', change_seed
    build_ref_hash_table(ref_kmer_dict, ref_file, k_size)
    # sort the node with the occurrence of the kmer
    node_lst = sorted(dbg_graph.nodes.values(), key=lambda x: x.occur, reverse=True)
    seeds = []
    kmer_mask = (1 << (k_size << 1)) - 1
    for one_node in node_lst:
        for seq_index in range(one_node.len - k_size + 1):
            # 从ref_bin_seq中获取指定长度的kmer
            bin_kmer = (one_node.value >> (seq_index << 1)) & kmer_mask
            if bin_kmer in ref_kmer_dict:
                seeds.append(one_node)
                break
        if len(seeds) >= seed_num:
            break
    if not seeds:
        seed_num = seed_num if len(node_lst) > seed_num else len(node_lst)
        seeds = node_lst[:seed_num]
    # Candidate seeds refer to nodes that frequently occur and whose sequences are also present within the reference sequences
    for seed in seeds:
        # traverse the graph from the seed
        path = dbg_graph.dfs_backtrack_traversal(seed, max_iterations, False) + [
            seed] + dbg_graph.dfs_backtrack_traversal(
            seed, max_iterations, True)
        # get the last node
        end_node = path.pop()
        # for node in path, the last k_size - 2 bases are the same with the first k_size - 2 bases of the next node
        # so cut the last k_size-2 bases and set the length of the node to len - k_size + 2
        path = [binseq2dna(node.value >> (2 * (k_size - 2)), node.len - k_size + 2) for node in path]
        seq = ''.join(path) + binseq2dna(end_node.value, end_node.len)
        if len(seq) > len(best_seq):
            best_seq = seq
            best_seed = seed

    # correct the orientation of the sequence
    kmer_sum = [0, 0]
    bin_seq_lst = [dnaseq2int(best_seq), dnaseq2int(dna_reverse_complement(best_seq))]
    for index, bin_seq_info in enumerate(bin_seq_lst):
        bin_seq, bin_seq_len = bin_seq_info
        for seq_index in range(bin_seq_len - k_size + 1):
            bin_kmer = (bin_seq >> (seq_index << 1)) & kmer_mask
            kmer_sum[index] += ref_kmer_dict[bin_kmer]
    if kmer_sum[0] < kmer_sum[1]:
        best_seq, best_seed = dna_reverse_complement(best_seq), dna_reverse_complement(best_seed)
    return best_seq, best_seed


# traverse the graph from the init_node and return the sequence
def assemble_seq_from_dbg(dbg_graph: DeBruijnGraph, max_iterations=1000):
    # select the node with the most occurrence as the seed
    k_size = dbg_graph.k_size
    best_seq, best_seed = '', ''
    seeds = sorted(dbg_graph.nodes.values(), key=lambda x: x.occur, reverse=True)[:10]
    for seed in seeds:
        # traverse the graph from the seed
        path = dbg_graph.dfs_backtrack_traversal(seed, max_iterations, False) + [
            seed] + dbg_graph.dfs_backtrack_traversal(
            seed, max_iterations, True)
        # get the last node
        end_node = path.pop()
        # for node in path, the last k_size - 2 bases are the same with the first k_size - 2 bases of the next node
        # so cut the last k_size-2 bases and set the length of the node to len - k_size + 2
        path = [binseq2dna(node.value >> (2 * (k_size - 2)), node.len - k_size + 2) for node in path]
        seq = ''.join(path) + binseq2dna(end_node.value, end_node.len)
        if len(seq) > len(best_seq):
            best_seq = seq
            best_seed = seed
    return best_seq


def assemble_reads_to_seq(gene_name: str, fa_1: str, fa_2: str, ref_file_path: str, assemble_kmer: int, output_dir: str,
                          kmer_limit: int = 4,
                          dynamic_kmer: bool = False):
    assemble_start_t = time.perf_counter()
    read_kmer_dict = defaultdict(int)
    # build the hash table based on the reads file
    logger.debug("Build the hash table based on the reads file!")
    build_read_hash_table(fa_1, fa_2, read_kmer_dict, assemble_kmer)
    # remove low occurrence kmers based on the kmer_limit
    read_kmer_dict = {kmer: count for kmer, count in read_kmer_dict.items() if count >= kmer_limit}
    dBG = DeBruijnGraph(k_size=assemble_kmer)
    # build the de bruijn graph
    logger.debug("Build the de bruijn graph and perform simplification and merging!")
    dBG.build_de_bruijn_graph(read_kmer_dict)
    dBG.simple_de_bruijn_graph()
    dBG.dbg_compress_nodes()
    dbg_build_t = time.perf_counter()
    logger.debug("The de bruijn graph has been built, used {:.2f}s!".format(dbg_build_t - assemble_start_t))
    contig, seed = assemble_seq_from_dbg_v2(dBG, ref_file_path)
    project_name = os.path.basename(output_dir.rstrip("/"))
    # contig = dna_reverse_complement(contig) if rev_comp else contig
    with open(os.path.join(output_dir, gene_name + ".fasta"), "wt") as outfile:
        seq_id = ">" + project_name + "_" + gene_name + "_k" + str(assemble_kmer) + '_l' + str(len(contig)) + "\n"
        outfile.write(seq_id + contig + "\n")
    assemble_end_t = time.perf_counter()
    logger.debug("The gene {} has been assembled, used {:.2f}s".format(gene_name, assemble_end_t - assemble_start_t))


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Easy353's Assembler
       zzhen@sculab ''')
    pars.add_argument('-i', dest='input_dir', help='Input a file(directory) with reads data', required=True)
    pars.add_argument("-r", dest="reference", type=str, help="Input a file(directory) with references", required=True)
    pars.add_argument("-o", dest="output_dir", type=str, help="Output directory.", required=False,
                      default="assemble_out")
    pars.add_argument("-ak", dest="assemble_kmer", type=int, help="Kmer setting for assembling reads. Default:31",
                      default=31)
    pars.add_argument("-at", dest="assemble_thread", type=int, help="Threads setting for assembling reads. Default:4",
                      default=4)
    pars.add_argument("-kmer_limit", dest="kmer_limit", type=int, help="Limit of kmer count. Default:3", default=2)
    pars.add_argument("-change_seed", dest="change_seed", type=int, help="Times of changing seed. Default:32",
                      default=32)
    pars.add_argument("-gdk", dest="get_dynamic_kmer", action="store_true",
                      help="Whether to use dynamic kmer as the kmer used for assembly.")
    pars.add_argument('-silent', dest='silent', action='store_true',
                      help='If using silent mode, there will not output any info!')
    pars.add_argument('-log', dest='log_file', type=str, help='The log file.', default='assemble.log')
    pars.add_argument('-paired', dest='paired_reads', action='store_true',
                      help='If using paired reads, please add this parameter.')
    args = pars.parse_args()
    # args.input_dir = "/Users/zzhen/Desktop/test/filter_out"
    # args.reference = "/Users/zzhen/Desktop/test/Apiaceae353"
    # args.output_dir = "/Users/zzhen/Desktop/test/assemble_out3"
    # args.assemble_kmer = 31
    # args.kmer_limit = 4
    # args.change_seed = 32
    # args.assemble_thread = 10
    # args.get_dynamic_kmer = False
    # args.log_file = 'assembly.log'
    # args.silent = False
    # args.paired_reads = True

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
    print_level = logging.ERROR if args.silent else logging.INFO
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
    t_assemble_start = time.perf_counter()

    ref_path_lst = get_file_lst(args.reference)
    # 读取reference文件
    if not ref_path_lst:
        logger.error('There were no reference files inputted.')
        sys.exit(-1)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    else:
        logger.warning('The output directory: {} exists, it will be deleted and remade!'.format(args.output_dir))
        shutil.rmtree(args.output_dir)
        os.makedirs(args.output_dir)
    logger.info('{} genes will be assembled with using {} threads.'.format(len(ref_path_lst), args.assemble_thread))

    if args.assemble_thread == 1:
        # when filter_thread is 1
        for ref_file_name in ref_path_lst:
            gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
            fa_1 = os.path.join(args.input_dir, gene_name + "_R1.fasta") if args.paired_reads else os.path.join(
                args.input_dir, gene_name + ".fasta")
            fa_2 = os.path.join(args.input_dir, gene_name + "_R2.fasta") if args.paired_reads else fa_1
            if not os.path.isfile(fa_1) or not os.path.isfile(fa_2):
                logger.error('The read files of {} were not found!'.format(gene_name))
                continue
            assemble_reads_to_seq(gene_name, fa_1, fa_2, ref_file_name, args.assemble_kmer, args.output_dir,
                                  args.kmer_limit, args.get_dynamic_kmer)
    else:
        with ProcessPoolExecutor(max_workers=min(args.assemble_thread, len(ref_path_lst))) as executor:
            tasks = []
            for ref_file_name in ref_path_lst:
                gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
                fa_1 = os.path.join(args.input_dir, gene_name + "_R1.fasta") if args.paired_reads else os.path.join(
                    args.input_dir, gene_name + ".fasta")
                fa_2 = os.path.join(args.input_dir, gene_name + "_R2.fasta") if args.paired_reads else fa_1
                if not os.path.isfile(fa_1) or not os.path.isfile(fa_2):
                    logger.error('The read files of {} were not found!'.format(gene_name))
                    continue
                tasks.append(
                    executor.submit(assemble_reads_to_seq, gene_name, fa_1, fa_2, ref_file_name, args.assemble_kmer,
                                    args.output_dir, args.kmer_limit, args.get_dynamic_kmer))
            task_count = 0
            for future in as_completed(tasks):
                task_count += 1
                if task_count % 10 == 0:
                    logger.info("{} / {} genes have been assembled!".format(task_count, len(ref_path_lst)))
                if task_count == len(ref_path_lst):
                    logger.info("All genes have been assembled!".format(task_count, len(ref_path_lst)))

    t_assemble_end = time.perf_counter()
    logger.info('Total time used for assembly: {:.2f} s'.format(t_assemble_end - t_assemble_start))
