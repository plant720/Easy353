# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/5/23 21:51
# @Author     : zzhen
# @File       : process_combined_seq.py
# @Software   : PyCharm
# @Description: 用来对于合并后的序列进行处理，包括：多序列比对MSA、序列修剪trim，构建基因树和物种树build gene trees and the species tree
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.

import os
import sys
import shutil
import logging
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger('EasyPiper')


# return a dict that stores the gene name and the path of the gene file
# 返回一个文件夹下所有以.fa .fasta等后缀结尾的文件的文件名和路径
def get_file_dict(in_dir: str) -> dict:
    file_path_dict = {}
    ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
    if not os.path.isdir(in_dir):
        logger.error("The path {} is not a directory!".format(in_dir))
    for i in os.listdir(in_dir):
        if os.path.isfile(os.path.join(in_dir, i)) and i.endswith(ext):
            file_path_dict[os.path.basename(i).split(".")[0]] = os.path.join(in_dir, i)
    return file_path_dict


# 进行多序列比对
def do_alignment(combined_file: str, aligned_file: str):
    cmd = 'muscle3 -in {} -out {} -quiet'.format(combined_file, aligned_file)
    subprocess.call(cmd, shell=True)
    return combined_file


def trim_alignment(aligned_file: str, trimmed_file: str):
    cmd = 'trimal -in {} -out {} -automated1'.format(aligned_file, trimmed_file)
    subprocess.call(cmd, shell=True)
    return trimmed_file


# build gene tree with raxmlHPC
def build_tree_raxml(trimmed_file: str, file_ext: str, raxml_result_dir: str):
    if not os.path.isdir(raxml_result_dir):
        os.makedirs(raxml_result_dir)
    # -n 文件名后缀 -w 输出文件夹 -o 外类群 -s 输入文件 -f bootstrap -m model -N bootstrap次数
    cmd = 'raxmlHPC -s {} -m GTRGAMMA -p 12345 -n {} -f a -x 12345 -N 100 -w {}'.format(trimmed_file, file_ext,
                                                                                        raxml_result_dir)
    subprocess.call(cmd, shell=True)
    # 返回raxmlHPC的结果文件
    return raxml_result_dir + '/RAxML_bipartitions.' + file_ext


# 将基因树文件合并到一个文件中
def combine_trees(tree_dir: str, tree_file: str):
    # 将所有的结果文件合并到一个文件中
    with open(tree_file, 'w') as f:
        for i in os.listdir(tree_dir):
            with open(os.path.join(tree_dir, i), 'r') as f1:
                f.writelines(f1.readlines())


# build species tree with astral
def build_species_tree(in_file: str, out_file: str):
    cmd = 'astral -i {} -o {}'.format(in_file, out_file)
    subprocess.call(cmd, shell=True)
    # 返回raxmlHPC的结果文件
    return out_file


# 检测是否存在命令
def check_command(cmd: str):
    return shutil.which(cmd)


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='EasyPiper zzhen@sculab ')
    # # 允许传入多个文件夹
    pars.add_argument('-i', metavar='<str>', type=str, help='The combined sequences dir', required=True)
    pars.add_argument('-o', metavar='<str>', type=str, help='The output dir', required=False)
    pars.add_argument('-t', metavar='<int>', type=int, default=10, help='Multi-threading number: default 10')
    pars.add_argument('-alignment', action='store_true', help='Whether to align the combine files with muscle')
    pars.add_argument('-trim', action='store_true', help='Whether to trim the alignment files with trimAl')
    pars.add_argument('-tree', action='store_true',
                      help='Whether to use multi-species coalescent-model to build the phylogenetic trees with RAxML and ASTRAL')
    args = pars.parse_args()

    # set the logger
    logger = logging.getLogger('EasyPiper')
    logger.setLevel(logging.DEBUG)

    # 建立一个filehandler来把日志记录在文件里，级别为debug以上
    fh = logging.FileHandler('EasyPiper.log', mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(lineno)s : %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S'))
    # 将相应的handler添加在logger对象中
    logger.addHandler(fh)
    # 建立一个stream-handler来把日志打在CMD窗口上，级别为error以上
    ch = logging.StreamHandler(sys.stdout)
    # 当使用silent mode时，只输出error日志
    print_level = logging.INFO
    ch.setLevel(print_level)
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)

    output_dir = args.o
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # 比对文件夹
    alignment_dir = os.path.join(output_dir, 'alignment_out')
    # 修剪文件夹
    trim_dir = os.path.join(output_dir, 'trim_out')
    # raxml_result
    raxml_dir = os.path.join(output_dir, 'raxmlHPC_out')
    # gene_tree 文件
    gene_tree_dir = os.path.join(output_dir, 'gene_tree')
    # astral 文件夹
    astral_dir = os.path.join(output_dir, 'astral_out')

    logger.info("=" * 80)
    logger.info("EasyPiper: A script to process the recovered combined genes (You can modify it for your own purpose)")
    logger.info("Author: zzhen")
    logger.info("The script can help you to align, trim, build gene trees and species trees")
    logger.info("The script is based on muscle, trimAl, raxmlHPC and astral. Please make sure they are available")
    logger.info(
        "Notice: This script serves as a straightforward workflow aimed at enhancing everyone's comprehension of the relevant processing steps. However, during actual analysis, it may be necessary to employ more intricate and nuanced approaches.")
    logger.info("=" * 80)

    if args.alignment:
        logger.info("1. Alignment start".center(80, '*'))
        if check_command('muscle3') is None:
            logger.error('The muscle3 for alignment is not available')
            exit(1)
        logger.info('The command for alignment is [muscle3 -in {} -out {} -quiet]'.format('<combine_file>',
                                                                                          '<alignment_file>'))
        if not os.path.isdir(alignment_dir):
            os.makedirs(alignment_dir)
        combine_dir_file_dict = get_file_dict(args.i)
        results = []
        with ThreadPoolExecutor(max_workers=args.t) as executor:
            futures = [
                executor.submit(do_alignment, combine_file, os.path.join(alignment_dir, gen_name + '.alignment.fasta'))
                for gen_name, combine_file in combine_dir_file_dict.items()]
            for future in as_completed(futures):
                results.append(future.result())
                if len(results) % 20 == 0:
                    logger.info('{} files have been aligned'.format(len(results)))
        logger.info('Alignment done!'.center(80, '*'))
    if args.trim:
        logger.info("2. Trim start".center(80, '*'))
        if not os.path.isdir(trim_dir):
            os.makedirs(trim_dir)
        # 检测是否存在trimal
        if check_command('trimal') is None:
            logger.error('The trimal for trimming is not available')
            exit(1)
        logger.info(
            "The command for trimming is [trimal -in {} -out {} -automated1]".format('<alignment_file>', '<trim_file>'))
        alignment_dir_file_dict = get_file_dict(alignment_dir)
        results = []
        with ThreadPoolExecutor(max_workers=args.t) as executor:
            futures = [
                executor.submit(trim_alignment, alignment_file, os.path.join(trim_dir, gen_name + '.trim.fasta'))
                for gen_name, alignment_file in alignment_dir_file_dict.items()]
            for future in as_completed(futures):
                results.append(future.result())
                if len(results) % 20 == 0:
                    logger.info('{} files have been trimmed'.format(len(results)))
        logger.info('Trim done!'.center(80, '*'))

    if args.tree:
        # 并联法 使用raxml and astral
        logger.info("3. Tree start".center(80, '*'))
        logger.info("This step includes building gene trees with RAxML and the species tree with ASTRAL")
        if check_command('raxmlHPC') is None:
            logger.error("The RAxML for building gene trees is not available")
            exit(1)
        logger.info(
            "The command for building gene trees is [raxmlHPC -s {} -m GTRGAMMA -p 12345 -n {} -f a -x 12345 -N 100 -w {}]".format(
                '<trim_file>', '<gen_name>', '<raxml_dir>'))
        if check_command('astral') is None:
            logger.error('The astral for building the species tree is not available')
            exit(1)
        logger.info("The command for building the species tree is [astral -i {} -o {}]".format('<gene_tree_file>',
                                                                                               '<species_tree_file>'))

        if not os.path.isdir(raxml_dir):
            os.makedirs(raxml_dir)
        if not os.path.isdir(gene_tree_dir):
            os.makedirs(gene_tree_dir)
        if not os.path.isdir(astral_dir):
            os.makedirs(astral_dir)
        trim_dir_file_dict = get_file_dict(trim_dir)
        results = []
        with ThreadPoolExecutor(max_workers=args.t) as executor:
            futures = [
                executor.submit(build_tree_raxml, trim_file, gen_name + '.tre',
                                os.path.join(os.getcwd(), raxml_dir, gen_name))
                for gen_name, trim_file in trim_dir_file_dict.items()]
            for future in as_completed(futures):
                results.append(future.result())
                if len(results) % 20 == 0:
                    logger.info('{} / {} gene trees have been built'.format(len(results), len(trim_dir_file_dict)))
        logger.info('Build gene trees done!'.center(80, '*'))
        # 将单个物种的基因树移动到一个文件夹中
        for result in results:
            shutil.copyfile(result, os.path.join(gene_tree_dir, result.split('/')[-1]))
        # 将基因树文件进行合并到一个文件夹中
        astral_in_file = os.path.join(astral_dir, 'combine.tre')
        combine_trees(gene_tree_dir, astral_in_file)
        # astral 输出文件
        astral_out_file = os.path.join(astral_dir, 'astral.tre')
        # astral 命令
        build_species_tree(astral_in_file, astral_out_file)
        logger.info('Build species tree done!'.center(80, '*'))
    shutil.copy("EasyPiper.log", os.path.join(output_dir, 'EasyPiper.log'))
