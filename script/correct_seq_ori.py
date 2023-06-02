# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/6/1 14:56
# @Author     : zzhen
# @File       : correct_seq_ori.py
# @Software   : PyCharm
# @Description: 用于校正序列的方向，使得序列的方向与参考序列的方向一致
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.


import os
import shutil
import argparse
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


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
                print(correct_seq_file)
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


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='EasyPiper zzhen@sculab ')
    # required parameters
    pars.add_argument('-i', metavar='<str>', type=str, help='The combined sequences dir', required=True)
    pars.add_argument('-r', metavar='<str>', type=str, help='The reference sequences dir', required=True)
    pars.add_argument('-o', metavar='<str>', type=str, help='The output dir', required=False, default="corrected_seqs")
    pars.add_argument("-t", metavar="<int>", type=int, default=10, help="Multi-threading number. Default: 10")
    args = pars.parse_args()

    multiprocess_correct_seq_ori(args.r, args.i, args.o, args.t)
