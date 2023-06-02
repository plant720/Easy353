# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/5/22 22:53
# @Author     : zzhen
# @File       : add_353species.py
# @Software   : PyCharm
# @Description: 添加treeoflife.kew.org官网上其他物种的353序列参与构建系统发育树
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.

import os
import shutil
import argparse
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


# generate gene reference from species 353 files
def generate_gene_ref_from_species(species_file_dir: str, out_dir: str):
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    species_file_lst = [os.path.join(species_file_dir, _file_) for _file_ in os.listdir(species_file_dir) if
                        _file_.endswith(".fasta")]
    for species_file in species_file_lst:
        species_name = os.path.basename(species_file).split(".")[-3]
        for record in SeqIO.parse(species_file, "fasta"):
            with open(os.path.join(out_dir, record.id + ".fasta"), "a") as out_file:
                out_file.write(">{}\n".format(species_name))
                out_file.write(str(record.seq) + "\n")


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


# return the reverse complementary sequence of a DNA sequence.
def dna_reverse_complement(_dna_: str) -> str:
    # use str.maketrans to build a translation table
    return _dna_.upper().translate(str.maketrans("ACGTU", "TGCAA", "RYMKSWHBVDN\n"))[::-1]


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


# build a hash table based on the reference file which is a fasta file
def build_ref_hash_table(ref_kmer_dict: defaultdict, ref_file_path: str, k_size: int):
    for record in SeqIO.parse(ref_file_path, "fasta"):
        ref_seq = "".join(filter(str.isalpha, record.seq.upper()))
        # convert the DNA sequence to a binary sequence
        for seq_index in range(len(ref_seq) - k_size + 1):
            # 从ref_bin_seq中获取指定长度的kmer
            kmer = ref_seq[seq_index:seq_index + k_size]
            if kmer in ref_kmer_dict:
                ref_kmer_dict[kmer] += 1
            else:
                ref_kmer_dict[kmer] = 1


def correct_seq_orientation(ref_file: str, ori_seqs_file: str, correct_seq_file: str):
    ref_kmer_dict = defaultdict(int)
    k_size = 21
    build_ref_hash_table(ref_kmer_dict, ref_file, k_size)
    ori_seq_info = []
    with open(correct_seq_file, "w") as out_file:
        for record in SeqIO.parse(ori_seqs_file, "fasta"):
            ori_seq = "".join(filter(str.isalpha, record.seq.upper()))
            forward_kmer = sum(
                [ref_kmer_dict[ori_seq[seq_index:seq_index + k_size]] for seq_index in
                 range(len(ori_seq) - k_size + 1)])
            reverse_kmer = sum(
                [ref_kmer_dict[dna_reverse_complement(ori_seq)[seq_index:seq_index + k_size]] for seq_index in
                 range(len(ori_seq) - k_size + 1)])
            out_file.write(">{}\n".format(record.id))
            if forward_kmer < reverse_kmer:
                out_file.write(str(dna_reverse_complement(ori_seq)) + "\n")
            else:
                out_file.write(ori_seq + "\n")
            ori_seq_info.append((record.id, forward_kmer, reverse_kmer))


# multiprocess
def multiprocess_correct_seq_ori(ref_seq_dir: str, ori_seqs_dir: str, correct_seqs_dir: str, thread_num: int = 8):
    if os.path.isdir(correct_seqs_dir):
        shutil.rmtree(correct_seqs_dir)
    if not os.path.isdir(correct_seqs_dir):
        os.makedirs(correct_seqs_dir)
    with ThreadPoolExecutor(max_workers=thread_num) as executor:
        tasks = []
        for gene_file_name in os.listdir(ori_seqs_dir):
            if not os.path.isfile(os.path.join(ref_seq_dir, gene_file_name)):
                continue
            tasks.append(executor.submit(correct_seq_orientation, os.path.join(ref_seq_dir, gene_file_name),
                                         os.path.join(ori_seqs_dir, gene_file_name),
                                         os.path.join(correct_seqs_dir, gene_file_name)))
        task_count = 0
        for future in as_completed(tasks):
            task_count += 1
            if task_count % 10 == 0:
                print("{} / {} files have had their sequences orientation corrected.".format(task_count, len(tasks)))
            if task_count == len(tasks):
                print("All {} files have had their sequences orientation corrected.".format(task_count))


# combine the add reference and the original reference
def add_more_seq_to_ori(ori_seq_file: str, add_seq_file: str, out_file: str):
    with open(out_file, "w") as out:
        for record in SeqIO.parse(ori_seq_file, "fasta"):
            out.write(">{}\n{}\n".format(record.id, record.seq))
        for record in SeqIO.parse(add_seq_file, "fasta"):
            out.write(">{}\n{}\n".format(record.id, record.seq))


def multiprocess_add_more_seq_to_ori(ori_seq_dir: str, add_seq_dir: str, out_seq_dir: str, thread_num: int = 8):
    if os.path.isdir(out_seq_dir):
        shutil.rmtree(out_seq_dir)
    if not os.path.isdir(out_seq_dir):
        os.makedirs(out_seq_dir)
    with ThreadPoolExecutor(max_workers=8) as executor:
        tasks = []
        for gene_file_name in os.listdir(add_seq_dir):
            if not os.path.isfile(os.path.join(ori_seq_dir, gene_file_name)):
                continue
            tasks.append(executor.submit(add_more_seq_to_ori, os.path.join(ori_seq_dir, gene_file_name),
                                         os.path.join(add_seq_dir, gene_file_name),
                                         os.path.join(out_seq_dir, gene_file_name)))
        task_count = 0
        for future in as_completed(tasks):
            task_count += 1
            if task_count % 10 == 0:
                print("{} / {} files have been combined.".format(task_count, len(tasks)))
            if task_count == len(tasks):
                print("All {} files have been combined".format(task_count))


if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='EasyPiper zzhen@sculab ')
    # required parameters
    pars.add_argument('-i', metavar='<str>', type=str,
                      help='The directory that stores the 353species files from treeoflife.kew.org', required=True)
    pars.add_argument('-c', metavar='<str>', type=str, help='The combined 353 seqs that recovered by Easy353',
                      required=True)
    pars.add_argument('-o', metavar='<str>', type=str, help='The output dir', required=False,
                      default="add_353species_out")
    pars.add_argument("-t", metavar="<int>", type=int, default=10, help="Multi-threading number. Default: 10")
    args = pars.parse_args()

    # 从treeoflife.kew.org下载的不同物种的353序列文件所在的文件夹
    species_file_dir = args.i
    # Easy353拼接并且合并后的序列文件
    combined_recovered_seqs = args.c
    # thread_num
    thread_num = args.t
    # 最终输出的文件夹
    out_dir = args.o
    if os.path.isdir(out_dir):
        print("The output dir {} already exists,it will be removed and remade!".format(out_dir))
        shutil.rmtree(out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    # 根据以物种存储的353序列的文件，生成一个文件夹，里面包含以基因作为分类的353序列
    add_353gene_dir = os.path.join(out_dir, "add_353genes")
    # 根据以基因作为分类的353序列，校正Easy353拼接产生的基因序列的方向
    correct_combined_seq = os.path.join(out_dir, "correct_combined_seqs")
    # 最终的合并后的序列
    final_combined_seq = os.path.join(out_dir, "final_combined_seqs")

    # 根据下载的多个物种的353基因序列，生成一个文件夹，里面包含每个基因的序列
    generate_gene_ref_from_species(species_file_dir, add_353gene_dir)
    # 根据下载的标准的基因序列，作为参考，用来校正Easy353拼接产生的基因序列的方向
    multiprocess_correct_seq_ori(ref_seq_dir=add_353gene_dir, ori_seqs_dir=combined_recovered_seqs,
                                 correct_seqs_dir=correct_combined_seq, thread_num=thread_num)
    # 将需要合并的序列合并到原始序列中
    multiprocess_add_more_seq_to_ori(ori_seq_dir=correct_combined_seq, add_seq_dir=add_353gene_dir, out_seq_dir=final_combined_seq,
                                     thread_num=thread_num)
