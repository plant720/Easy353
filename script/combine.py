#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/23 11:30 上午
# @Author     : zzhen
# @File       : combine.py
# @Software   : PyCharm
# @Description: 用于将多个物种的结果进行合并到一个文件中
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import os
import argparse
from collections import defaultdict


# 用于获取每个基因的统计信息：即该基因在哪几个物种中存在
# segment: contig or scaffold
def statistics_gene_info(dir_list: list, segment: str = "contig") -> dict:
    # 统计每一个物种恢复的基因数量
    species_info = defaultdict(list)
    for _dir_ in dir_list:
        genes_list = [i.replace("." + segment + ".fasta", "") for i in os.listdir(os.path.join(_dir_, segment))]
        species_info[_dir_] = genes_list
    # 统计恢复了某一个基因的物种
    genes = []
    genes_info = defaultdict(list)
    for _, value in species_info.items():
        genes.extend(value)
    genes = list(set(genes))
    for key, value in species_info.items():
        for gene in genes:
            if gene in value:
                genes_info[gene].append(key)
    return genes_info


def combine_contig_file(dir_list: list, out_dir: str, species_num, percent: float = 0.8,
                        segment: str = "contig") -> None:
    if not os.path.isdir(os.path.join(out_dir, segment)):
        os.makedirs(os.path.join(out_dir, segment))
    genes_info = statistics_gene_info(dir_list, segment)
    # 将每个物种的contig文件合并到一个文件中
    for gene, species_list in genes_info.items():
        # 如果该基因恢复了某几个物种，则将该基因的所有物种的contig文件合并到一个文件中
        if len(species_list) >= species_num * percent:
            with open(os.path.join(out_dir, segment, gene + "." + segment + ".fasta"), "a") as f:
                for species in species_list:
                    with open(os.path.join(species, segment, gene + "." + segment + ".fasta")) as f1:
                        for line in f1:
                            f.write(line)


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''MinerScript 
        zzhen@sculab ''')
    # 允许传入多个文件夹
    pars.add_argument('-f', metavar='<str>', type=str, help='''the directory of results''', required=True,
                      nargs="+")
    pars.add_argument('-o', metavar='<str>', type=str, help='''the output dir''', required=False,
                      default="combine")
    pars.add_argument("-p", metavar='<float>', type=float, help='''the percent of species''', required=False,
                      default=0.8)
    pars.add_argument("-scaffold", action="store_true", help="whether to combine scaffold")
    args = pars.parse_args()

    # 一个存放多个文件夹的列表
    results_dir = args.f
    scaffold = args.scaffold
    output_dir = args.o
    percent = args.p
    combine_contig_file(results_dir, output_dir, species_num=len(results_dir), percent=percent, segment="contig")
    if scaffold:
        combine_contig_file(results_dir, output_dir, species_num=len(results_dir), percent=percent, segment="scaffold")
