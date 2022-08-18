#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/7/31 2:04 下午
# @Author     : zzhen
# @File       : combine.py.py
# @Software   : PyCharm
# @Description: 将Easy353做出的结果文件 进行合并 并进行比对 并进行trim
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import collections
import glob
import os
import argparse
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
from Bio import SeqIO


# 返回文件夹下的所有基因文件
def get_file_dict(_dir_path_: str) -> dict:
    _file_path_dict_ = {}
    if os.path.isdir(_dir_path_):
        for _ in glob.glob(os.path.join(_dir_path_, '*')):
            ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
            if os.path.isfile(_) and _.endswith(ext):
                _file_path_dict_[_.split('/')[-1].split(".")[0]] = _
    else:
        print('Error: The input path is not a directory!')
    return _file_path_dict_


# 用于记录每个物种的统计信息：即该物种拼接出的基因名称
def statistics_species(dir_list: list):
    ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
    gene_info = defaultdict(list)
    for _dir_ in dir_list:
        genes_list = [i.replace(".target.fasta", "") for i in os.listdir(os.path.join(_dir_, "target_genes")) if
                      i.endswith(ext)]
        gene_info[_dir_] = genes_list
    return gene_info


# 用于记录每个基因的统计信息：即该基因在哪几个物种中存在
def statistics_gene(dir_list: list) -> dict:
    species_info = statistics_species(dir_list)
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


# 根据所有物种的共有的基因，将每个物种的target文件合并到一个文件中
def combine_target_gene(dir_list: list, combine_dir: str, species_num, percent: float = 1.0) -> None:
    genes_info = statistics_gene(dir_list)
    # 将每个物种的contig文件合并到一个文件中
    for gene, species_list in genes_info.items():
        # 如果该基因恢复了某几个物种，则将该基因的所有物种的contig文件合并到一个文件中
        if len(species_list) >= species_num * percent:
            with open(os.path.join(combine_dir, gene + ".target.fasta"), "a") as f:
                for species in species_list:
                    with open(os.path.join(species, "target_genes", gene + ".target.fasta")) as f1:
                        f.write(f1.readline().split("_")[0] + "\n")
                        f.writelines(f1.readlines())


# 进行多序列比对
def do_alignment(combined_file: str, aligned_file: str):
    cmd = "muscle3 -in {} -out {} -quiet".format(combined_file, aligned_file)
    subprocess.call(cmd, shell=True)
    return combined_file


def trim_alignment(aligned_file: str, trimmed_file: str):
    cmd = "trimal -in {} -out {} -automated1".format(aligned_file, trimmed_file)
    subprocess.call(cmd, shell=True)
    return trimmed_file


# build gene tree with raxmlHPC
def build_tree_raxml(trimmed_file: str, file_ext: str, raxml_result_dir: str):
    if not os.path.isdir(raxml_result_dir):
        os.makedirs(raxml_result_dir)
    # -n 文件名后缀 -w 输出文件夹 -o 外类群 -s 输入文件 -f bootstrap -m model -N bootstrap次数
    cmd = "raxmlHPC -s {} -m GTRGAMMA -p 12345 -n {} -f a -x 12345 -N 1000 -w {}".format(trimmed_file, file_ext,
                                                                                         raxml_result_dir)
    subprocess.call(cmd, shell=True)
    # 返回raxmlHPC的结果文件
    return raxml_result_dir + "/RAxML_bipartitions." + file_ext


# 将基因树文件合并到一个文件中
def combine_trees(tree_dir: str, tree_file: str):
    # 将所有的结果文件合并到一个文件中
    with open(tree_file, "w") as f:
        for i in os.listdir(tree_dir):
            with open(os.path.join(tree_dir, i), "r") as f1:
                f.writelines(f1.readlines())


# build species tree with astral
def build_species_tree(in_file: str, out_file: str):
    cmd = "astral -i {} -o {}".format(in_file, out_file)
    subprocess.call(cmd, shell=True)
    # 返回raxmlHPC的结果文件
    return out_file


# 将多个alignment file合并成一个文件
def combine_alignment(alignment_file_dir: str, combine_file: str):
    ext = (".fasta", ".fa", ".fas", ".fsa", ".faa")
    alignment_files = [os.path.join(alignment_file_dir, i) for i in os.listdir(alignment_file_dir) if i.endswith(ext)]
    # read each seq from fasta file
    seq_dict = collections.defaultdict(str)
    for record in SeqIO.parse(alignment_files[0], "fasta"):
        seq_dict[record.id] = str(record.seq)
    for alignment_file in alignment_files[1:]:
        for record in SeqIO.parse(alignment_file, "fasta"):
            seq_dict[record.id] += str(record.seq)
    # write seq to fasta file
    with open(combine_file, "w") as f:
        for key, value in seq_dict.items():
            f.write(">" + key + "\n" + value + "\n")


# 使用mrbayes 进行建树
def build_species_tree_mrbayes(mrbayes_in_file: str, mrbayes_out_file: str):
    pass


# 检测是否存在命令
def check_command(cmd: str):
    return shutil.which(cmd)


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''MinerScript 
        zzhen@sculab ''')
    # 允许传入多个文件夹
    pars.add_argument('-i', metavar='<str>', type=str, help='''the directory of results''', required=True,
                      nargs="+")
    pars.add_argument('-o', metavar='<str>', type=str, help='''the output dir''', required=False,
                      default="combine")
    pars.add_argument("-p", metavar='<float>', type=float, help='''the percent of species''', required=False,
                      default=1.0)
    pars.add_argument("-t", metavar="<int>", type=int, default=10, help="How many threads are used in alignment")
    pars.add_argument("-alignment", action="store_true", help="Whether to align the combine files")
    pars.add_argument("-trim", action="store_true", help="Whether to trim the alignment files")
    pars.add_argument("-tree", action="store_true", help="Whether to build the phylogenetic tree with RAxML")
    pars.add_argument("-concatenate", action="store_true",
                      help="Whether to use concatenated-model to build the phylogenetic trees")
    pars.add_argument("-coalescent", action="store_true",
                      help="Whether to use multi-species coalescent-model to build the phylogenetic trees")

    args = pars.parse_args()
    # 一个存放多个文件夹的列表
    results_dir = args.i
    output_dir = args.o
    species_percent = args.p
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "combine_stats.csv")

    out_stream = open(output_file, "a")
    # 统计每个物种拼接出的基因名称
    out_stream.write("Genes recovered from each species:\n")
    for key, value in statistics_species(results_dir).items():
        out_stream.write(key + "," + ",".join(value) + "\n")
    # 统计每个基因被哪几个物种拼接出
    out_stream.write("*************************************************************\n")
    out_stream.write("*************************************************************\n")
    out_stream.write("Each gene is recovered by which species:\n")
    for key, value in statistics_gene(results_dir).items():
        out_stream.write(key + "," + ",".join(value) + "\n")
    out_stream.close()

    # 合并文件夹
    combine_dir = os.path.join(output_dir, "combine")
    if not os.path.isdir(combine_dir):
        os.makedirs(combine_dir)
    # 比对文件夹
    alignment_dir = os.path.join(output_dir, "alignment")
    if not os.path.isdir(alignment_dir):
        os.makedirs(alignment_dir)
    # 修剪文件夹
    trim_dir = os.path.join(output_dir, "trim")
    if not os.path.isdir(trim_dir):
        os.makedirs(trim_dir)
    # raxml_result
    raxml_dir = os.path.join(output_dir, "raxmlHPC")
    if not os.path.isdir(raxml_dir):
        os.makedirs(raxml_dir)
    # gene_tree 文件
    gene_tree_dir = os.path.join(output_dir, "gene_tree")
    if not os.path.isdir(gene_tree_dir):
        os.makedirs(gene_tree_dir)
    # astral 文件夹
    astral_dir = os.path.join(output_dir, "astral")
    if not os.path.isdir(astral_dir):
        os.makedirs(astral_dir)
    # mrbayes 文件夹
    mrbayes_dir = os.path.join(output_dir, "mrbayes")
    if not os.path.isdir(mrbayes_dir):
        os.makedirs(mrbayes_dir)

    combine_target_gene(results_dir, combine_dir, species_num=len(results_dir), percent=species_percent)
    print("Combine target genes done!")

    if args.alignment:
        # 检测是否存在muscle3
        if check_command("muscle3") is None:
            print("The muscle3 for alignment is not available")
            exit(1)
        combine_dir_file_dict = get_file_dict(combine_dir)
        results = []
        with ThreadPoolExecutor(max_workers=args.t) as executor:
            futures = [
                executor.submit(do_alignment, combine_file, os.path.join(alignment_dir, gen_name + ".alignment.fasta"))
                for gen_name, combine_file in combine_dir_file_dict.items()]
            for future in as_completed(futures):
                results.append(future.result())
                if len(results) % 20 == 0:
                    print("{} files have been aligned".format(len(results)))
        print("Alignment done!")
    if args.trim:
        # 检测是否存在trimal
        if check_command("trimal") is None:
            print("The trimal for trim is not available")
            exit(1)
        alignment_dir_file_dict = get_file_dict(alignment_dir)
        results = []
        with ThreadPoolExecutor(max_workers=args.t) as executor:
            futures = [
                executor.submit(trim_alignment, alignment_file, os.path.join(trim_dir, gen_name + ".trim.fasta"))
                for gen_name, alignment_file in alignment_dir_file_dict.items()]
            for future in as_completed(futures):
                results.append(future.result())
                if len(results) % 20 == 0:
                    print("{} files have been trimmed".format(len(results)))
        print("Trim done!")

    if args.tree:
        # 并联法 使用raxml and astral
        if args.coalescent:
            # 检测是否存在RAxML 和 astral
            if check_command("raxmlHPC") is None or check_command("astral") is None:
                print("The RAxML and astral are not available")
                exit(1)
            trim_dir_file_dict = get_file_dict(trim_dir)
            results = []
            with ThreadPoolExecutor(max_workers=args.t) as executor:
                futures = [
                    executor.submit(build_tree_raxml, trim_file, gen_name + ".tre",
                                    os.path.join(os.getcwd(), raxml_dir, gen_name))
                    for gen_name, trim_file in trim_dir_file_dict.items()]
                for future in as_completed(futures):
                    results.append(future.result())
                    if len(results) % 20 == 0:
                        print("{} trees have been built".format(len(results)))
            print("Build gene trees done!")
            # 将单个物种的基因树移动到一个文件夹中
            for result in results:
                shutil.copyfile(result, os.path.join(gene_tree_dir, result.split("/")[-1]))
            # 将基因树文件进行合并到一个文件夹中
            astral_in_file = os.path.join(astral_dir, "combine.tre")
            combine_trees(gene_tree_dir, astral_in_file)
            # astral 输出文件
            astral_out_file = os.path.join(astral_dir, "astral.tre")
            # astral 命令
            build_species_tree(astral_in_file, astral_out_file)
            print("Build species tree done!")
        # 串联法 使用mrbayes
        if args.concatenate:
            # 检测是否存在mrbayes
            if check_command("mb") is None:
                print("The mrbayes is not available")
                exit(1)
            trim_dir_file_dict = get_file_dict(trim_dir)
            # 将多个基因修剪好的alignment合并到一个文件中
            mrbayes_in_file = os.path.join(mrbayes_dir, "combine.alignment.fasta")
            combine_alignment(trim_dir, mrbayes_in_file)
            mrbayes_out_file = os.path.join(mrbayes_dir, "mrbayes.tre")
            # build_species_tree_mrbayes(mrbayes_in_file, mrbayes_out_file)
