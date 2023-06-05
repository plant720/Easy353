#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/3/18 3:44 下午
# @Author     : zzhen
# @File       : easy353.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import logging
import argparse
import multiprocessing
import os
import shutil
import sys
import platform
import time
from collections import defaultdict
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed

from Easy353Lib.filter import make_ref_kmer_dict, reads_bin_filter, get_file_lst
from Easy353Lib.assemble import assemble_reads_to_seq

if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage='%(prog)s [options]',
                                   description="Easy353 zzhen@sculab")
    pars.add_argument('-1', dest='fq_file_1', type=str, help='Specify the input sequencing file(s) -1 (.fq/.gz).',
                      required=False)
    pars.add_argument('-2', dest='fq_file_2', type=str, help='Specify the input sequencing file(s) -2 (.fq/.gz).',
                      required=False)
    pars.add_argument('-r', dest='reference', type=str, help='Specify an input reference file or directory.')
    pars.add_argument('-o', dest='out_dir', type=str, help='Designate the output directory.', default='easy353_out')
    pars.add_argument('-fk', dest='filter_kmer', type=int, help='Set the kmer for filtering reads. Default:21',
                      default=21)
    pars.add_argument('-s', dest='step_length', type=int,
                      help='Determine the length of the sliding window on the reads. Default:1', default=1)
    pars.add_argument('-ft', dest='filter_thread', type=int,
                      help='Set the number of threads for filtering reads. Default:4', default=4)
    pars.add_argument("-ak", dest="assemble_kmer", type=int, help="Set the kmer for assembling reads. Default:31",
                      default=31)
    pars.add_argument("-at", dest="assemble_thread", type=int, help="Set the number of threads for assembling reads. Default:8",
                      default=8)
    pars.add_argument("-kmer_limit", dest="kmer_limit", type=int, help="Set a limit of kmer count. Default:4", default=4)
    pars.add_argument('-fpr', dest='filter_pair_read', action='store_true',
                      help='Enable fpr mode to get more filtered reads. Tradeoff between obtaining more reads (yes) or saving time (no).')
    pars.add_argument("-change_seed", dest="change_seed", type=int, help="Specify the number of seed changes. Default:32",
                      default=32)
    pars.add_argument("-gdk", dest="get_dynamic_kmer", action="store_true",
                      help="Choose whether to use dynamic kmer as the kmer used for assembly.")
    pars.add_argument('-log', dest='log_file', type=str, help='Specify the log file.', default='easy353.log')
    pars.add_argument('-silent', dest='silent', action='store_true',
                      help='Enable silent mode. When active, the program will not output any information.')
    args = pars.parse_args()

    # args.fq_file_1 = '/Volumes/zzhen/data/work_data/Easy353/empirical_data/tran_apiaceae/Anthriscus_sylvestris.1.fq.gz'
    # args.fq_file_2 = '/Volumes/zzhen/data/work_data/Easy353/empirical_data/tran_apiaceae/Anthriscus_sylvestris.2.fq.gz'
    # args.reference = '/Users/zzhen/Desktop/test/Apiaceae353'
    # args.out_dir = '/Users/zzhen/Desktop/test/Apiaceae353_out'
    # args.filter_kmer = 21
    # args.step_length = 1
    # args.filter_thread = 8
    # args.filter_pair_read = False
    # args.silent = False
    # args.assemble_kmer = 31
    # args.kmer_limit = 4
    # args.change_seed = 32
    # args.assemble_thread = 10
    # args.get_dynamic_kmer = False
    # args.log_file = 'easy353.log'
    # args.silent = False

    # set the logger
    logger = logging.getLogger("easy353")
    logger.setLevel(logging.DEBUG)

    # 建立一个stream-handler来把日志打在CMD窗口上，级别为error以上
    ch = logging.StreamHandler(sys.stdout)
    # 当使用silent mode时，只输出error日志
    print_level = logging.WARNING if args.silent else logging.INFO
    ch.setLevel(print_level)
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)

    # 先检查文件夹是否都存在
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    else:
        # 输出文件夹不允许存在
        logger.warning('The output directory: {} exists, it will be deleted and remade!'.format(args.out_dir))
        shutil.rmtree(args.out_dir)
        os.makedirs(args.out_dir)

    # 建立一个filehandler来把日志记录在文件里，级别为debug以上
    fh = logging.FileHandler(os.path.join(args.out_dir, args.log_file))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(lineno)s : %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S'))
    # 将相应的handler添加在logger对象中
    logger.addHandler(fh)

    #  set the start method of multiprocessing
    if platform.system() in ('Linux', 'Darwin'):
        # fork方式 子进程可以共享父进程的资源
        multiprocessing.set_start_method('fork')
    else:
        # spawn方式 子进程会copy一份父进程的资源，导致内存占用较大
        multiprocessing.set_start_method('spawn')
    logger.info("="*80)
    logger.info("Easy353: A tool to get Angiosperms353 genes for phylogenomic research")
    logger.info("Version: 2.0.1")
    logger.info("Author: zzhen")
    logger.info("fq_file_1         : {}".format(args.fq_file_1))
    logger.info("fq_file_2         : {}".format(args.fq_file_2))
    logger.info("reference         : {}".format(args.reference))
    logger.info("out_dir           : {}".format(args.out_dir))
    logger.info("filter_kmer       : {}".format(args.filter_kmer))
    logger.info("filter_thread     : {}".format(args.filter_thread))
    logger.info("step_length       : {}".format(args.step_length))
    logger.info("assemble_kmer     : {}".format(args.assemble_kmer))
    logger.info("assemble_thread   : {}".format(args.assemble_thread))
    logger.info("kmer_limit        : {}".format(args.kmer_limit))
    logger.info("filter_pair_read  : {}".format(args.filter_pair_read))
    logger.info("change_seed       : {}".format(args.change_seed))
    logger.info("get_dynamic_kmer  : {}".format(args.get_dynamic_kmer))
    logger.info("log_file          : {}".format(args.log_file))
    logger.info("silent            : {}".format(args.silent))
    logger.info("=" * 80)

    # 检查参数是否合理
    if not os.path.isfile(args.fq_file_1) and not os.path.isfile(args.fq_file_2):
        logger.error("The fq_file {} and {} do not exist!".format(args.fq_file_1, args.fq_file_2))
        sys.exit(1)
    if not os.path.isdir(args.reference):
        logger.error("The reference directory {} does not exist!".format(args.reference))
        sys.exit(1)
    if args.filter_thread < 1:
        logger.error("The filter_thread must be greater than 0!")
        sys.exit(1)
    if args.assemble_thread < 1:
        logger.error("The assemble_thread must be greater than 0!")
        sys.exit(1)
    if args.step_length < 1:
        logger.error("The step_length must be greater than 0!")
        sys.exit(1)
    if args.kmer_limit < 1:
        logger.error("The kmer_limit must be greater than 0!")
        sys.exit(1)
    if args.change_seed < 1:
        logger.error("The change_seed must be greater than 0!")
        sys.exit(1)

    filter_out = os.path.join(args.out_dir, 'filter_out')
    if not os.path.isdir(filter_out):
        os.makedirs(filter_out)
    assemble_out = os.path.join(args.out_dir, 'assemble_out')
    if not os.path.isdir(assemble_out):
        os.makedirs(assemble_out)
    logger.info(" Easy353 is running ".center(80, "="))
    logger.info(" 1. Build reference hash table ".center(80, "="))

    t_hash_start = time.perf_counter()
    ref_kmer_dict = defaultdict(int)
    make_ref_kmer_dict(ref_kmer_dict, args.reference, args.filter_kmer)
    logger.info('Hash dictionary has been made.')
    t_hash_end = time.perf_counter()
    logger.info('Time used for building hash table: {:.2f} s'.format(t_hash_end - t_hash_start))
    logger.info('The memory usage of the current process is: {:.2f} GB'.
                format(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024))
    logger.info(" 2. Read filtering ".center(80, "="))
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
                         filter_out, 0, args.filter_thread, args.filter_pair_read)
    else:
        process_list = []
        for p_id in range(args.filter_thread):
            p = multiprocessing.Process(target=reads_bin_filter,
                                        args=(args.reference, ref_kmer_dict, args.filter_kmer, args.step_length,
                                              args.fq_file_1, args.fq_file_2,
                                              filter_out, p_id, args.filter_thread, args.filter_pair_read,))
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
            with open(os.path.join(filter_out, prefix + ".fasta"), 'wt') as merge_file:
                for file_name in [os.path.join(filter_out, prefix + "." + str(i + 1) + ".fasta") for i in
                                  range(args.filter_thread)]:
                    with open(file_name, 'rt') as infile:
                        shutil.copyfileobj(infile, merge_file)
                    os.remove(file_name)

    t_filter_end = time.perf_counter()
    logger.info('Time used for filter: {:.2f} s'.format(t_filter_end - t_hash_end))
    logger.info('Total time used for filter: {:.2f} s'.format(t_filter_end - t_hash_start))

    t_assemble_start = time.perf_counter()
    logger.info(' 3. Read Assembly '.center(80, "="))
    logger.info('{} genes will be assembled with using {} threads.'.format(len(ref_path_lst), args.assemble_thread))
    if args.assemble_thread == 1:
        # when filter_thread is 1
        for ref_file_name in ref_path_lst:
            gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
            fa_1 = os.path.join(filter_out, gene_name + "_R1.fasta") if paired_reads else os.path.join(
                filter_out, gene_name + ".fasta")
            fa_2 = os.path.join(filter_out, gene_name + "_R2.fasta") if paired_reads else fa_1
            assemble_reads_to_seq(gene_name, fa_1, fa_2, ref_file_name, args.assemble_kmer, assemble_out,
                                  args.kmer_limit, args.get_dynamic_kmer)
    else:
        with ProcessPoolExecutor(max_workers=min(args.assemble_thread, len(ref_path_lst))) as executor:
            tasks = []
            for ref_file_name in ref_path_lst:
                gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
                fa_1 = os.path.join(filter_out, gene_name + "_R1.fasta") if paired_reads else os.path.join(
                    filter_out, gene_name + ".fasta")
                fa_2 = os.path.join(filter_out, gene_name + "_R2.fasta") if paired_reads else fa_1
                tasks.append(executor.submit(assemble_reads_to_seq, gene_name, fa_1, fa_2, ref_file_name,
                                             args.assemble_kmer,
                                             assemble_out, args.kmer_limit, args.get_dynamic_kmer))
            count = 0
            for future in as_completed(tasks):
                count += 1
                if count % 10 == 0:
                    logger.info("{} / {} genes have been assembled!".format(count, len(ref_path_lst)))
                if count == len(ref_path_lst):
                    logger.info("All genes have been assembled!".format(count, len(ref_path_lst)))
    t_assemble_end = time.perf_counter()
    logger.info('Total time used for assembly: {:.2f} s'.format(t_assemble_end - t_assemble_start))
    logger.info('All jobs have been done!')
    logger.info('Total time used: {:.2f} s'.format(t_assemble_end - t_hash_start))
