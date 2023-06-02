# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/5/20 13:33
# @Author     : zzhen
# @File       : multi_sample_gene_merge.py
# @Software   : PyCharm
# @Description: 用于将多个样本恢复的Angiosperms353基因合并到一个文件中
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.


import os
import csv
import shutil
import argparse
from collections import defaultdict


# return a dict that stores the gene name and the path of the gene file
# 返回一个文件夹下所有以.fa .fasta等后缀结尾的文件的文件名和路径
def get_file_dict(in_dir: str) -> dict:
    file_path_dict = {}
    ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
    if not os.path.isdir(in_dir):
        raise ValueError("The path {} is not a directory!".format(in_dir))
    for i in os.listdir(in_dir):
        if os.path.isfile(os.path.join(in_dir, i)) and i.endswith(ext):
            file_path_dict[os.path.basename(i).split(".")[0]] = os.path.join(in_dir, i)
    return file_path_dict


def dict_to_csv(data: dict, filename: str):
    # 获取字典的所有键作为列名
    columns = list(data.keys())

    # 获取字典中所有值的行名
    rows = set()
    for values in data.values():
        rows.update(values.keys())

    # 将字典输出为CSV文件
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[''] + columns)
        writer.writeheader()
        for row in rows:
            row_data = {column: data[column].get(row, '') for column in columns}
            writer.writerow({'': row, **row_data})


# 记录每个样本恢复的基因名称
def sample_gene_recovery_statistics(samples_result_lst: list, assemble_dir: str = "assemble_out"):
    ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
    sample_gene_recovery_info_dict = defaultdict(list)
    for sample_name in samples_result_lst:
        genes_list = [i.replace(".fasta", "") for i in os.listdir(os.path.join(sample_name, assemble_dir)) if
                      i.endswith(ext)]
        # 处理名称，只保留样本名，去除绝对路径
        sample_name = os.path.basename(os.path.abspath(sample_name))
        sample_gene_recovery_info_dict[sample_name] = genes_list
    return sample_gene_recovery_info_dict


# 用于记录每个基因的统计信息：即该基因在哪几个物种中存在
def gene_distribution(samples_result_lst: list) -> dict:
    species_info_dict = sample_gene_recovery_statistics(samples_result_lst)
    genes = []
    genes_distribution_dict = defaultdict(list)
    for _, value in species_info_dict.items():
        genes.extend(value)
    genes = list(set(genes))
    for key, value in species_info_dict.items():
        for gene in genes:
            if gene in value:
                genes_distribution_dict[gene].append(key)
    return genes_distribution_dict


# 根据所有物种的共有的基因，将每个物种的target文件合并到一个文件中
def multi_sample_gene_merge(samples_result_lst: list, combine_dir: str, species_num, percent: float = 1.0) -> None:
    genes_distribution_dict = gene_distribution(samples_result_lst)
    # 将每个物种的contig文件合并到一个文件中
    for gene, species_list in genes_distribution_dict.items():
        # 如果该基因恢复了某几个物种，则将该基因的所有物种的contig文件合并到一个文件中
        if len(species_list) >= species_num * percent:
            with open(os.path.join(combine_dir, gene + ".fasta"), "a") as f:
                for species in species_list:
                    with open(os.path.join(species, "assemble_out", gene + ".fasta")) as f1:
                        next(f1)
                        f.write(">" + species + "\n")
                        f.writelines(f1.readlines())


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='EasyPiper zzhen@sculab ')
    pars.add_argument('-i', metavar='<str>', type=str,
                      help='The result directories of some samples after running Easy353', required=True, nargs="+")
    pars.add_argument('-o', metavar='<str>', type=str, help='The output dir', required=False, default="combine")
    pars.add_argument("-p", metavar="<float>", type=float, default=1.0,
                      help="The percent of species that a gene is recovered,default=1.0")
    pars.add_argument("-t", metavar="<int>", type=int, default=10, help="Multi-threading number,default=10")
    args = pars.parse_args()

    samples_result_lst = args.i
    output_dir = args.o
    thread_num = args.t
    percent = args.p
    if os.path.isdir(output_dir):
        print("The output dir {} is already exist, it will be removed and remade!".format(output_dir))
        shutil.rmtree(output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # set the log file
    log_file = os.path.join(output_dir, "combine.log")
    csv_file = os.path.join(output_dir, "combine.csv")

    with open(log_file, "w") as f:
        f.write("Combine the recovered genes from different samples!\n")
        samples_result_lst = [os.path.basename(os.path.abspath(i)) for i in samples_result_lst]
        genes_distribution_dict = gene_distribution(samples_result_lst)
        f.write("The result directories of samples: {}\n".format(",".join(samples_result_lst)))
        f.write("The recovered gene names: {}\n".format(",".join(genes_distribution_dict.keys())))
        f.write("The genes that are recovered by more than {}% species will be combined!\n".format(percent * 100))
        f.write("The recovered genes are combined in the directory: {}\n".format(os.path.abspath(output_dir)))
        sample_and_gene_full_info_dict = {}
        for gene_name, species_list in genes_distribution_dict.items():
            tmp = {i: (i in species_list) for i in samples_result_lst}
            sample_and_gene_full_info_dict[gene_name] = tmp
        f.write("The recovered genes and their distribution in each sample will be written in {}:\n".format(csv_file))
        dict_to_csv(sample_and_gene_full_info_dict, csv_file)
        f.write("Each gene is recovered by which sample:\n")
        for key, value in genes_distribution_dict.items():
            f.write(key + "," + ",".join(value) + "\n")
        f.write("=="*30+"\n")
        f.write("The recovered genes from each sample:\n")
        for key, value in sample_gene_recovery_statistics(samples_result_lst).items():
            f.write(key + "," + ",".join(value) + "\n")
        # 合并文件夹
        combine_dir = os.path.join(output_dir, "combine")
        if not os.path.isdir(combine_dir):
            os.makedirs(combine_dir)
        multi_sample_gene_merge(samples_result_lst, combine_dir, species_num=len(samples_result_lst), percent=percent)
