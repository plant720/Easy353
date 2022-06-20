#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/3/18 3:44 下午
# @Author     : zzhen
# @File       : easy353.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import argparse
from src import assemble, filter


# 参数初始化
def args_init():
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage="%(prog)s [options]",
                                   description="Easy353 zzhen@sculab")
    pars.add_argument("-1", dest="fq_file_1", type=str, nargs="+",
                      help="Input file(s) with forward paired-end reads (*.fq/.gz/.tar.gz).", required=False)
    pars.add_argument("-2", dest="fq_file_2", type=str, nargs="+",
                      help="Input file(s) with reverse paired-end reads (*.fq/.gz/.tar.gz).", required=False)
    pars.add_argument("-u", dest="unpaired_fq_file", type=str,
                      help="Input file(s) with unpaired (single-end) reads.", required=False, nargs="+")
    pars.add_argument("-r", dest="reference", type=str, help="Input a file(directory) with references.", required=True)
    pars.add_argument("-o", dest="output_dir", type=str, help="Output directory.", required=False,
                      default="easy353_output")
    pars.add_argument("-k1", dest="filter_kmer", type=int, help="Kmer setting for filtering reads. Default:31",
                      default=31)
    pars.add_argument("-k2", dest="assemble_kmer", type=int, help="Kmer setting for assembling reads. Default:41",
                      default=41)
    pars.add_argument("-s", dest="step_length", type=int,
                      help="Step length of the sliding window on the reads. Default:1", default=1)
    pars.add_argument("-t1", dest="filter_thread", type=int,
                      help="Threads setting for filtering reads. Default:4", default=4)
    pars.add_argument("-t2", dest="assemble_thread", type=int,
                      help="Threads setting for assembling reads. Default:4", default=4)
    pars.add_argument("-kmer_limit", dest="kmer_limit", type=int, help="Limit of kmer count. Default:2", default=2)
    pars.add_argument("-f", dest="function_mode", type=int, help="0:all,1:filter,2:assemble. Default:0", default=0)
    pars.add_argument("-min", dest="minimum_length_ratio", type=float,
                      help="The minimum ratio of contig length to reference average length. Default:1.0", default=1.0)
    pars.add_argument("-max", dest="maximum_length_ratio", type=float,
                      help="The maximum ratio of contig length to reference average length. Default:2.0", default=2.0)
    pars.add_argument("-change_seed", dest="change_seed", type=int, help="Times of changing seed. Default:32",
                      default=32)
    pars.add_argument("-scaffold", dest="generate_scaffold", action="store_true",
                      help="Whether to generate scaffolds or not.")
    pars.add_argument("-fast", dest="fast", action="store_true", help="Whether to use fast mode.")
    # 默认限制使用的参考序列的数量
    pars.add_argument("-reference_limit", dest="reference_limit", action="store_false",
                      help="Whether to limit the numbers of reference sequences or not.")
    args = pars.parse_args()
    return args


if __name__ == "__main__":
    args = args_init()
    # python easy353.py -1 /Users/zzhen/Desktop/ara.1.fq -2 /Users/zzhen/Desktop/ara.2.fq -r /Users/zzhen/Desktop/new -o /Users/zzhen/Desktop/test_a -k1 31 -k2 41 -s 1 -t1 1 -t2 1
    fastq_files = tuple()
    _paired_reads_ = False
    if args.unpaired_fq_file:
        fastq_files = (args.unpaired_fq_file,)
    if args.fq_file_1 and args.fq_file_2:
        fastq_files = (args.fq_file_1, args.fq_file_2)
        _paired_reads_ = True
    if args.function_mode == 0 or args.function_mode == 1:
        filter.filter_flow(_read_data_tuple_=fastq_files, _out_dir_=args.output_dir, _reference_path_=args.reference,
                           _kmer_size_=args.filter_kmer, _step_size_=args.step_length,
                           _ref_reverse_complement_=args.fast, _paired_reads_=_paired_reads_,
                           _thread_for_filter_=args.filter_thread)
    if args.function_mode == 0 or args.function_mode == 2:
        assemble.assemble_flow(_input_read_path_=args.output_dir, _out_dir_=args.output_dir, _ref_path_=args.reference,
                               _assemble_kmer_size_=args.assemble_kmer, _assemble_thread_=args.assemble_thread,
                               _ref_reverse_complement_=True, _pos_=True,
                               _change_seed_=args.change_seed, _kmer_limit_count_=args.kmer_limit,
                               _min_percent_length_=args.minimum_length_ratio,
                               _max_percent_length_=args.maximum_length_ratio,
                               _iteration_=1000, _write_scaffold_=args.generate_scaffold)
    # out_dir = "/Users/zzhen/Desktop/result"
    # reads_dir = "reads"
    # re_assemble_dir = "re-assemble"
    # if not os.path.isdir(os.path.join(out_dir, reads_dir)):
    #     os.makedirs(os.path.join(out_dir, reads_dir))
    # if not os.path.isdir(os.path.join(out_dir, re_assemble_dir)):
    #     os.makedirs(os.path.join(out_dir, re_assemble_dir))
    # # 已经成功组装的基因名
    # assembled_gene = [i.split(".")[0] for i in os.listdir(os.path.join(out_dir, "contig")) if i.endswith(".fasta")]
    # for i in os.listdir(out_dir):
    #     path = os.path.join(out_dir, i)
    #     if os.path.isfile(path) and path.endswith(".fasta"):
    #         if i.split(".")[0] in assembled_gene:
    #             shutil.move(path, os.path.join(out_dir, reads_dir, i))
    #         else:
    #             shutil.move(path, os.path.join(out_dir, re_assemble_dir, i))
    # print("Done!")
