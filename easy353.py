#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/3/18 3:44 下午
# @Author     : zzhen
# @File       : easy353.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import argparse
import multiprocessing
import os
import sys
import platform
from Easy353Lib import assemble, filter


def main():
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage="%(prog)s [options]",
                                   description="Easy353 zzhen@sculab")
    pars.add_argument("-1", dest="fq_file_1", type=str, nargs="+",
                      help="Input file(s) with forward paired-end reads (*.fq/.gz/.tar.gz).")
    pars.add_argument("-2", dest="fq_file_2", type=str, nargs="+",
                      help="Input file(s) with reverse paired-end reads (*.fq/.gz/.tar.gz).")
    pars.add_argument("-u", dest="unpaired_fq_file", type=str,
                      help="Input file(s) with unpaired (single-end) reads.", nargs="+")
    pars.add_argument("-r", dest="reference", type=str, help="Input a file(directory) with references.")
    pars.add_argument("-o", dest="output_dir", type=str, help="Output directory.", default="easy353_output")
    pars.add_argument("-k1", dest="filter_kmer", type=int, help="Kmer setting for filtering reads. Default:29",
                      default=29)
    pars.add_argument("-k2", dest="assemble_kmer", type=int, help="Kmer setting for assembling reads. Default:41",
                      default=41)
    pars.add_argument("-s", dest="step_length", type=int,
                      help="Step length of the sliding window on the reads. Default:1", default=1)
    pars.add_argument("-t1", dest="filter_thread", type=int,
                      help="Threads setting for filtering reads. Default:1", default=1)
    pars.add_argument("-t2", dest="assemble_thread", type=int,
                      help="Threads setting for assembling reads. Default:4", default=4)
    pars.add_argument("-kmer_limit", dest="kmer_limit", type=int, help="Limit of kmer count. Default:8", default=8)
    pars.add_argument("-f", dest="function_mode", type=int, help="0:all,1:filter,2:assemble. Default:0", default=0)
    pars.add_argument("-min", dest="minimum_length_ratio", type=float,
                      help="The minimum ratio of contig length to reference average length. Default:1.0", default=1.0)
    pars.add_argument("-max", dest="maximum_length_ratio", type=float,
                      help="The maximum ratio of contig length to reference average length. Default:2.0", default=2.0)
    pars.add_argument("-change_seed", dest="change_seed", type=int, help="Times of changing seed. Default:32",
                      default=32)
    # 默认限制使用的参考序列的数量
    pars.add_argument("-reference_number", dest="reference_number", type=int,
                      help="The number of the reference sequences used to build hash table. Default:all", default=None)
    pars.add_argument("-fast", dest="fast", action="store_true", help="Whether to use fast mode.")
    args = pars.parse_args()
    if len(sys.argv) <= 2:
        pars.print_help()
        sys.exit(0)
    fastq_files = tuple()
    _paired_reads_ = False
    if args.unpaired_fq_file:
        fastq_files = (args.unpaired_fq_file,)
    if args.fq_file_1 and args.fq_file_2:
        fastq_files = (args.fq_file_1, args.fq_file_2)
        _paired_reads_ = True
    filtered_read_dir = os.path.join(args.output_dir, "filtered_reads")
    if not os.path.isdir(filtered_read_dir):
        os.makedirs(filtered_read_dir)
    if args.function_mode == 0 or args.function_mode == 1:
        filter.filter_flow(_read_data_tuple_=fastq_files, _out_dir_=filtered_read_dir, _reference_path_=args.reference,
                           _kmer_size_=args.filter_kmer, _step_size_=args.step_length,
                           _ref_reverse_complement_=args.fast, _paired_reads_=_paired_reads_,
                           _thread_for_filter_=args.filter_thread, _ref_number_=args.reference_number)
    if args.function_mode == 0 or args.function_mode == 2:
        assemble.assemble_flow(_input_read_path_=filtered_read_dir, _out_dir_=args.output_dir,
                               _ref_path_=args.reference,
                               _assemble_kmer_size_=args.assemble_kmer, _assemble_thread_=args.assemble_thread,
                               _ref_reverse_complement_=True, _pos_=True,
                               _change_seed_=args.change_seed, _kmer_limit_count_=args.kmer_limit,
                               _min_percent_length_=args.minimum_length_ratio,
                               _max_percent_length_=args.maximum_length_ratio,
                               _iteration_=1000)


if __name__ == "__main__":
    if platform.system() in ("Linux", "Darwin"):
        multiprocessing.set_start_method('fork')
    else:
        multiprocessing.set_start_method('spawn')
    main()

