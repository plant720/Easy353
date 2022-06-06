#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/5/5 5:09 下午
# @Author     : zzhen
# @File       : statistics.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import os
import argparse
from collections import defaultdict


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Unsupported value encountered.')


# 用于统计所有物种共有的基因
def statistics(dir_list: list, _scaffold_: bool):
    ori_elements = [i.replace(".contig.fasta", "") for i in os.listdir(os.path.join(dir_list[0], "contig"))]
    if _scaffold_:
        ori_elements.extend([i.replace(".scaffold.fasta", "")
                             for i in os.listdir(os.path.join(dir_list[0], "scaffold"))])
        ori_elements = list(set(ori_elements))
    for _dir_ in dir_list[1:]:
        elements = [i.replace(".contig.fasta", "") for i in os.listdir(os.path.join(_dir_, "contig"))]
        if _scaffold_:
            elements.extend([i.replace(".scaffold.fasta", "")
                             for i in os.listdir(os.path.join(_dir_, "scaffold"))])
            elements = list(set(elements))
        ori_elements = list(set(ori_elements).intersection(set(elements)))
    return sorted(ori_elements)


# 用于记录每个物种的统计信息：即该物种拼接出的基因名称
def statistics_detail(dir_list: list, _scaffold_: bool):
    gene_info = defaultdict(list)
    for _dir_ in dir_list:
        genes_list = [i.replace(".contig.fasta", "") for i in os.listdir(os.path.join(_dir_, "contig"))]
        if _scaffold_:
            genes_list.extend([i.replace(".scaffold.fasta", "") for i in os.listdir(os.path.join(_dir_, "scaffold"))])
            genes_list = list(set(genes_list))
        gene_info[_dir_] = genes_list
    return gene_info


# 用于记录每个基因的统计信息：即该基因在哪几个物种中存在
def statistics_gene(dir_list: list, _scaffold_: bool) -> dict:
    species_info = statistics_detail(dir_list, _scaffold_)
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


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''MinerScript 
        zzhen@sculab ''')
    # 允许传入多个文件夹
    pars.add_argument('-f', metavar='<str>', type=str, help='''the directory of results''', required=True,
                      nargs="+")
    pars.add_argument('-o', metavar='<str>', type=str, help='''the output dir''', required=False,
                      default="combine")
    pars.add_argument("-scaffold", action="store_true", help="whether to combine scaffold")
    args = pars.parse_args()

    # 一个存放多个文件夹的列表
    results_dir = args.f
    scaffold = args.scaffold
    output_file = args.o
    out_stream = open(output_file, "a")

    # 统计每个物种拼接出的基因名称
    out_stream.write("the genes that each specie has in the results:\n")
    for key, value in statistics_detail(results_dir, scaffold).items():
        out_stream.write(key + "," + ",".join(value) + "\n")

    # 统计所有物种都拼接出的基因名称
    out_stream.write("*************************************************************\n")
    out_stream.write("*************************************************************\n")
    out_stream.write("the shared gene for all species in the results:\n")
    out_stream.write(",".join(statistics(results_dir, scaffold)) + "\n")

    # 统计每个基因被哪几个物种拼接出
    out_stream.write("*************************************************************\n")
    out_stream.write("*************************************************************\n")
    out_stream.write("the gene for which species in the results:\n")
    for key, value in statistics_gene(results_dir, scaffold).items():
        out_stream.write(key + "," + ",".join(value) + "\n")

    out_stream.close()
