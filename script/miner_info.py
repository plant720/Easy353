#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/21 8:45 上午
# @Author     : zzhen
# @File       : miner_info.py.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import os
import sys
from Bio import SeqIO

if __name__ == '__main__':
    data_dir = sys.argv[1]
    species_dir = [os.path.join(data_dir, x) for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, x))]
    species_log_csv = open(os.path.join(data_dir, 'species_log.csv'), 'a+')
    log_csv = open(os.path.join(data_dir, 'log.csv'), 'a+')
    for species in species_dir:
        if os.path.isdir(os.path.join(species, "contig")):
            contigs = [os.path.join(species, "contig", x) for x in os.listdir(os.path.join(species, "contig"))]
            # 计算过滤成功的条数
            species_log_csv.write(species + ',' + str(len(contigs)) + '\n')
            total_gene_len = 0
            for contig in contigs:
                record = SeqIO.read(contig, "fasta")
                gene_len = record.description.split("_")[-1]
                total_gene_len += int(gene_len)
                log_csv.write(contig + ',' + gene_len + '\n')
            species_log_csv.write(species + ',' + str(total_gene_len) + '\n')
    log_csv.close()
    species_log_csv.close()