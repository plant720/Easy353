#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/3/6 15:58
# @Author  : xiepulin
# @File    : my_filter.py
# @Software: PyCharm

import sys
import argparse
from Bio import SeqIO,Seq
import gzip
import os
import re
import time
from collections import defaultdict
import random
import  multiprocessing as mp
import  platform
import gc
import multiprocessing
from multiprocessing import Process,Pool
import shutil


'''
判断输入参考基因组是文件还是文件夹，0代表文件夹，1代表文件
'''
def file_or_directory(ref):
    if os.path.isdir(ref):
        files = os.listdir(ref)
        if len(files) == 0:
            print("the input  reference genome folder is empty")
            sys.exit()
        else:
            flag = 0
            # print("this is a dir")
    elif os.path.isfile(ref):
        size = os.path.getsize(ref)
        if size == 0:
            print("The input reference genome file is empty")
            sys.exit()
        else:
            flag = 1
            # print("this is a file")
    else:
        sys.exit()
    return flag


def get_absolute(path):
    if path==None:
        return None
    else:
        if os.path.isabs(path):
            return path
        else:
            path = os.path.abspath(path)
            return path

def get_file_physical_size(path):
    size=os.path.getsize(path)
    return size

def get_platform():
    current_os = platform.system().lower()   #linux windows darwin
    return  current_os


'''
返回参考序列列表   参考基因组可能是文件，也可能是文件夹
'''
def get_file_list(path):
    file_list=[]
    if os.path.isdir(path):
        files=get_files(path)
        for i in files:
            if os.path.getsize(i)==0:
                print("{} is empty".format(i))
                sys.exit()
            else:
                file_list.append(i)
    elif os.path.isfile(path):
        size = os.path.getsize(path)
        if size == 0:
            print("{}：is empty".format(path))
            sys.exit()
        else:
            file_list.append(path)
    else:
        sys.exit()
    return file_list


'''
获取目录下所有文件
'''
def get_files(ref):
    file_path_list = []
    for root, dirs, files in os.walk(ref):
        for file in files:
            file_path = os.path.join(root, file)
            file_path_list.append(file_path)
    return file_path_list

'''
文件是否真实存在
'''
def is_exist(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0

    elif os.path.isdir(file):
        files = get_files(file)
        if files==[]:
            flag=0
        else:
            flag = 1
            for i in files:
                if os.path.getsize(i) > 0:
                    continue
                else:
                    flag = 0
                    break
    else:
        flag=0

    return flag

'''
获取fasta文件的文件名 即基因名
'''
def get_basename(file):
    if is_exist(file):
        basename=os.path.basename(file)

        if ".fasta" in basename:
            basename = basename.split(".fasta")[0]
        elif ".fas" in basename:
            basename = basename.split(".fas")[0]
        elif ".fa" in basename:
            basename = basename.split(".fa")[0]
        else:
            basename=basename
        return basename


#获取反向互补
def reverse_complement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   return str("".join(complement.get(base, base) for base in reversed(seq)))

# 简化反向互补
def reverse_complement_limit(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

# 创建哈西字典,将reference存入字典,并记录位置信息
def gethashdict(reference, merSize, get_rc = False, pos = False,print_info=False):
    files = get_file_list(reference)
    if files == []:
        return 0
    gene_number = 0
    kmer_dict = defaultdict(list)
    for file in files:
        if os.path.isfile(file) == False:
            continue
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        file_name = get_basename(file)
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            gene_number += 1
            refseq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            for j in range(0, len(refseq)-merSize+1):
                temp_list, kmer = [], refseq[j:j+merSize]
                if kmer in kmer_dict:
                    temp_list = kmer_dict[kmer]
                    temp_list[0] += 1
                    if pos:
                        temp_list[1] += (j+1)/len(refseq)
                    if file_name not in temp_list:
                        temp_list.append(file_name)
                else:
                    temp_list = [1, (j+1)/len(refseq), file_name] if pos else [1, -1, file_name]
                    kmer_dict[sys.intern(kmer)] = temp_list
            if get_rc:
                try:
                    refseq =reverse_complement_limit(refseq)
                except:
                    print(file, refseq)
                for j in range(0, len(refseq)-merSize+1):
                    temp_list, kmer = [], refseq[j:j+merSize]
                    if kmer in kmer_dict:
                        temp_list = kmer_dict[kmer]
                        temp_list[0] += 1
                        if pos: temp_list[1] += (j+1)/len(refseq)
                        if file_name not in temp_list: temp_list.append(file_name)
                    else:
                        temp_list = [1, (j+1)/len(refseq), file_name] if pos else [1, -1, file_name]
                        kmer_dict[sys.intern(kmer)] = temp_list
            if print_info:
                print('Memory:', round(sys.getsizeof(kmer_dict)/1024/1024, 4), 'MB, Number of Genes:', gene_number, end='\r',flush=True)
        infile.close()
    return  kmer_dict




def filter_reads(_dict, _kmer, stepSize, readseq, get_rc = False):
    tmplist = []
    for j in range(0, len(readseq) - _kmer + 1 - stepSize, stepSize):
        kmer = readseq[j:j + _kmer]
        if kmer in _dict:
            tmplist.extend(_dict.get(kmer)[2:]) #剔除kmercount 和 位置信息

    if get_rc:   #data1的反向互补并不是data2 最正确的是过滤data1,过滤data2,取交集---------实际上根据data1过滤后的title定位data2也是可以的
        # temp_seq = Seq(readseq).reverse_complement()   #Bio.Seq(seq).reverse_complement()
        temp_seq = Seq.reverse_complement(readseq)       #Bio.Seq.reverse_complement(seq) 快
        for j in range(0, len(temp_seq) - _kmer + 1 - stepSize, stepSize):
            kmer = temp_seq[j:j + _kmer]
            if kmer in _dict:
                tmplist.extend(_dict.get(kmer)[2:])

    return set(tmplist)


def filter_paired_reads(_dict, _kmer,setpSize, file_1, file_2, out_dir, t_id, t_count, get_rc = False,data_size='all'):
    t1, t2, reads_count  = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    infile_2 = gzip.open(file_2, 'r') if file_2[-3:].lower() == ".gz" else open(file_2, 'rb')
    temp_rec1 = [infile_1.readline(),infile_1.readline(),infile_1.readline(),infile_1.readline()]
    temp_rec2 = [infile_2.readline(),infile_2.readline(),infile_2.readline(),infile_2.readline()]

    data_size= 1000000000000 if data_size.lower()=='all' else int(data_size)
    for _ in infile_1:
        reads_count += 1
        if reads_count% t_count == t_id:
            for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'), get_rc):   #如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
                with open(os.path.join(out_dir, file_name + ".fasta"), "a+") as outfile:
                    outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])
                    outfile.writelines(['>', temp_rec2[0].decode('utf-8')[1:], temp_rec2[1].decode('utf-8')])

        temp_rec1 = [_, infile_1.readline(),infile_1.readline(),infile_1.readline()]
        temp_rec2 = [infile_2.readline(),infile_2.readline(),infile_2.readline(),infile_2.readline()]


        if reads_count * t_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            # print('handled\t',reads_count * t_count //1000000, 'm reads, ',round(t2,2),'s/m reads', end='\r',sep="",flush=True)
            info = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count * t_count // 1000000, round(t2, 2))
            print(info, end='\r', flush=True)


        if reads_count * t_count >= data_size:
            break

    infile_1.close()
    infile_2.close()


def filter_single_reads(_dict, _kmer,setpSize, file_1,out_dir, t_id, t_count, get_rc = False,data_size='all'):
    t1, t2, reads_count = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    temp_rec1 = [infile_1.readline(), infile_1.readline(), infile_1.readline(), infile_1.readline()]

    data_size = 10000000000 if data_size.lower() == 'all' else int(data_size)
    for _ in infile_1:
        reads_count += 1
        if reads_count % t_count == t_id:
            for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'),
                                          get_rc):  # 如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
                with open(os.path.join(out_dir, file_name + ".fasta"), "a+") as outfile:
                    outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])

        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]

        if reads_count * t_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            print('handled\t', reads_count * t_count // 1000000, 'm reads, ', round(t2, 2), 's/m reads', end='\r',
                  sep="", flush=True)
        if reads_count * t_count >= data_size:
            break

    infile_1.close()



###################################
###################################
'''
re-filter
'''
###################################
###################################
'''
获得ref 平均长度
返回 {filename:length}
'''
def get_ref_info(reference,ref_length_dict):
    files = get_file_list(reference)
    for file in files:
        length_all = 0
        ref_count = 0
        file_name = get_basename(file)
        for rec in SeqIO.parse(file, "fasta"):
            length = len(rec.seq)
            length_all += length
            ref_count += 1
        average_length = int(length_all / ref_count)
        ref_length_dict[file_name] = average_length
    return ref_length_dict


'''
获取reads 平均长度 读取fastq文件
'''
def get_reads_length(file):
    infile = gzip.open(file, 'rt') if file[-3:].lower() == ".gz" else open(file, 'r')
    temp_rec = [infile.readline(), infile.readline(), infile.readline(), infile.readline()]

    reads_length = []
    line_number = 0
    line_limit = 100
    for _ in infile:
        if line_number < line_limit:
            reads_length.append(len(temp_rec[1]))
            line_number += 1
            temp_rec = [_, infile.readline(), infile.readline(),
                        infile.readline()]  # 每次读取的第一行 如果换成infile.readline(),相当于在此基础上又读了一行
        else:
            break
    infile.close()

    length = int(sum(reads_length) / len(reads_length))
    return length



'''
按照reads平均长度计算filtered_out 总bp数目，不影响结果但比详细计算bp数目快10倍
re-filter的结果为fasta格式 两行算一条read
'''
def get_reads_bp_number(filtered_out,reads_length,filtered_reads_whole_bp_dict):
    files = get_file_list(filtered_out)
    for file in files:
        file_name = get_basename(file)
        f=open(file,"rb")
        number=0
        for _ in f:
            number+=1
        number=int(number/2)
        whole_bp=reads_length*number
        filtered_reads_whole_bp_dict[file_name]=whole_bp
    return filtered_reads_whole_bp_dict


def re_filter_reads(_dict, _kmer, stepSize, readseq, get_rc=False):
    '''
    :param _dict:     ref中一个基因的hashdict
    :param _kmer:     wordsize
    :param stepSize:  step
    :param readseq:   read
    :param get_rc:    是否需要反向互补，filter中已经反向互补了，并写进一个文件。ref已经反向互补了，此时不需要反向互补了
    :return:
    '''
    tmplist = []
    for j in range(0, len(readseq) - _kmer + 1 - stepSize, stepSize):
        kmer = readseq[j:j + _kmer]
        if kmer in _dict:
            tmplist.extend(_dict.get(kmer)[2:])#剔除kmercount 和 位置信息

    if get_rc:
        # temp_seq = Seq(readseq).reverse_complement()   #Bio.Seq(seq).reverse_complement()
        temp_seq = Seq.reverse_complement(readseq)  # Bio.Seq.reverse_complement(seq) 快
        for j in range(0, len(temp_seq) - _kmer + 1 - stepSize, stepSize):
            kmer = temp_seq[j:j + _kmer]
            if kmer in _dict:
                tmplist.extend(_dict.get(kmer)[2:])

    return set(tmplist)


def do_re_filter_reads(parameter_information, used_dict, kmer, big_reads_path, re_filtered_reads_whole_bp_dict,
                       get_rc=False):
    '''
    :param parameter_information: 参数表
    :param used_dict: ref中一个基因的hashdict
    :param kmer: kmer
    :param big_reads_path: filtered_out超过500层或大于10M，big reads不发生变化
    :param re_filtered_reads_whole_bp_dict: 用于记录新生成的re-filter 过滤的reads bp总数
    :param get_rc: 是否需要反向互补
    :return:
    '''
    stepSize = parameter_information["step"]
    out_dir = parameter_information["out_dir"]

    infile = open(big_reads_path, "rb")
    temp_rec = [infile.readline(), infile.readline()]
    for _ in infile:
        for file_name in re_filter_reads(used_dict, kmer, stepSize, temp_rec[1].decode('utf-8'), get_rc):
            with open(out_dir + "/" + file_name + ".fasta", "a+") as f:
                re_filtered_reads_whole_bp_dict[file_name] += len(temp_rec[1].decode('utf-8'))
                f.writelines(['>', temp_rec[0].decode('utf-8')[1:], temp_rec[1].decode('utf-8')])
        temp_rec = [_, infile.readline()]

    infile.close()


def do_re_filter_loop(parameter_information, re_filtered_reads_whole_bp_dict, ref_length_dict, re_filter_dict,
                      thread_id, thread):
    '''
    :param parameter_information:
    :param re_filtered_reads_whole_bp_dict:
    :param ref_length_dict: 记录参考序列的平均长度，用于计算richness
    :param re_filter_dict: 记录需要重过滤基因 大于500层或大于10M
    :param thread_id:  标记线程id
    :param thread:     线程总数
    :return:
    '''

    current = 0
    for i in re_filter_dict:
        current += 1
        if current % thread == thread_id:
            big_reads_path = i["big_reads_path"]
            filtered_reads_path = i["filtered_reads_path"]
            filtered_reads_path_backup = i["filtered_reads_path_backup"]
            ref_path = i["ref_path"]
            gene_name = i["gene_name"]
            kmer = parameter_information["wordsize"]
            kmer_add = 2
            re_filter_number=0
            while True:
                kmer+=kmer_add
                if kmer>127:
                    break

                #备份数据 （每次开始优先备份数据，覆盖式的备份; 达到循环退出条件时候，删除备份）
                if re_filter_number==0:
                    shutil.copy(big_reads_path, filtered_reads_path_backup)
                else:
                    shutil.copy(filtered_reads_path, filtered_reads_path_backup)
                re_filter_number+=1

                used_dict = gethashdict(ref_path, kmer,get_rc=False,pos=False,print_info=False)
                with open(filtered_reads_path, "w") as f:  # 清空序列
                    pass
                re_filtered_reads_whole_bp_dict[gene_name] = 0  # 重置
                do_re_filter_reads(parameter_information, used_dict, kmer, big_reads_path,
                                   re_filtered_reads_whole_bp_dict,
                                   False)
                richness = int(re_filtered_reads_whole_bp_dict[gene_name] / ref_length_dict[gene_name])
                size = get_file_physical_size(filtered_reads_path)

                #相关参数
                sharp_cutoff_threshold=512*1024   #锐减阈值
                tolerate_threshold=10*1024*1024  #容忍阈值
                richness_level_1=512
                richness_level_2 = 1024
                richness_level_3 = 2048
                kmer_add_1=2
                kmer_add_2=4
                kmer_add_3=8


                #for test
                # sharp_cutoff_threshold = 200 * 1024  # 锐减阈值
                # tolerate_threshold = 250 * 1024  # 容忍阈值
                # richness_level_1 = 100
                # richness_level_2 = 1024
                # richness_level_3 = 2048
                # kmer_add_1 = 4
                # kmer_add_2 = 8
                # kmer_add_3 = 12

                #打印输出
                if size >=sharp_cutoff_threshold:
                    print("gene_name:{}  kmer:{} richness:{} size:{}KB".format(gene_name, kmer, richness,
                                                                               int(size / 1024)))

                #kmer太大，结果锐减且突破最低临界值 ---保留上一次的结果
                if size < sharp_cutoff_threshold:
                    os.remove(filtered_reads_path)
                    shutil.copy(filtered_reads_path_backup, filtered_reads_path)
                    os.remove(filtered_reads_path_backup)
                    break
                elif  size>=sharp_cutoff_threshold and size<=tolerate_threshold:
                    os.remove(filtered_reads_path_backup)
                    break
                else:
                    pass


                #确定Kmer增大幅度
                if richness <= richness_level_1 or kmer > 127 or size <= tolerate_threshold:
                    os.remove(filtered_reads_path_backup)
                    break
                elif richness >= richness_level_1 and richness < richness_level_2:
                    kmer_add = kmer_add_1
                    os.remove(filtered_reads_path_backup)
                elif richness >= richness_level_2 and richness <= richness_level_3:
                    kmer_add = kmer_add_2
                    os.remove(filtered_reads_path_backup)
                else:
                    kmer_add = kmer_add_3
                    os.remove(filtered_reads_path_backup)




                del used_dict
                gc.collect()  # 释放内存







def my_filter_main(args):
    # 初始化操作
    data1 = args.data1
    data2 = args.data2
    single = args.single
    thread = args.thread
    wordsize = args.wordsize
    out_dir = args.out
    step = args.step
    reference = args.reference
    data_size=args.data_size

    #读写瓶颈 硬盘
    if thread >= 4:
        limit_thread = 4
    else:
        limit_thread = thread



    parameter_information={"data1":data1,"data2":data2,"single":single,
                           "thread":thread,"wordsize":wordsize,"out_dir":out_dir,"step":step,
                           "reference":reference,"data_size":data_size
    }

    # print(parameter_information)

    print('======================== Filter =========================')
    hashdict = {}
    reads_length=100 #默认，后续会重新计算
    t1 = time.time()
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    hashdict = gethashdict(reference, wordsize,get_rc=True,pos=False,print_info=True)
    print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
    print("Hashdict has prepared")
    print('Memory of hashdict:', round(sys.getsizeof(hashdict) / 1000 / 1000, 4), 'MB')
    if sys.platform in ('darwin'):
        multiprocessing.set_start_method('fork')

    if data1 and data2:
        reads_length=get_reads_length(data1)
        thread_list = []
        for t in range(limit_thread):
            p = Process(target=filter_paired_reads,
                        args=(hashdict, wordsize,step,data1,data2,out_dir,t,limit_thread,False,data_size))
            thread_list.append(p)
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()

    if single:
        reads_length = get_reads_length(single)
        thread_list = []
        for t in range(limit_thread):
            p = Process(target=filter_single_reads,
                        args=(hashdict, wordsize, step, single, out_dir, t,limit_thread, False, data_size))
            thread_list.append(p)
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()

    t2=time.time()
    print(" "*50,flush=True,end='\r')  #之前的信息比较长，冲刷掉
    print("Time used: {}s".format(round((t2-t1),4)))

    ############################## re-filter #############################3
    print('======================== Re-Filter =========================')
    ref_length_dict = {}  # 参考序列的平均长度
    filtered_reads_whole_bp_dict = defaultdict(int)  # 第一次过滤的结果
    re_filtered_reads_whole_bp_dict = defaultdict(int)  ## 重过滤reads总数
    ref_length_dict = get_ref_info(reference, ref_length_dict)  # ref_length_dict  第一次过滤的参考信息
    filtered_reads_whole_bp_dict = get_reads_bp_number(out_dir, reads_length, filtered_reads_whole_bp_dict)
    re_filter_dict = []  # [{'gene_name': 'matK', 'richness': 627, 'size': 1418262, 'filtered_reads_path':path1,ref_path:path2}]
    if not os.path.isdir(os.path.join(out_dir, 'big_reads')):
        os.mkdir(os.path.join(out_dir, 'big_reads'))
    for i in zip(ref_length_dict.items(), filtered_reads_whole_bp_dict.items()):
        temp = {}
        gene_name = i[0][0]
        filtered_reads_path = os.path.join(out_dir, gene_name + ".fasta")
        filtered_reads_path_backup = os.path.join(out_dir, gene_name + "_backup" + ".fasta")
        if not is_exist(filtered_reads_path):
            continue
        if file_or_directory(reference):
            ref_path = reference
        else:
            ref_path = os.path.join(reference, gene_name + ".fasta")
        richness = int(i[1][1] / i[0][1])
        size = get_file_physical_size(filtered_reads_path)

        # if richness <= 100 or size <= 500000:  #for test
        #     continue
        if richness <= 512 or size <= 10 * 1024 * 1024:  # 低于512层或小于10MB 2**20
            continue
        else:
            big_reads_path = os.path.join(out_dir, "big_reads", gene_name + ".fasta")
            shutil.move(filtered_reads_path, big_reads_path)
            temp["gene_name"] = gene_name
            temp["richness"] = richness
            temp["size"] = size
            temp["filtered_reads_path"] = filtered_reads_path
            temp["filtered_reads_path_backup"] = filtered_reads_path_backup
            temp["ref_path"] = ref_path
            temp["big_reads_path"] = big_reads_path
            re_filter_dict.append(temp)
    if len(re_filter_dict)>0:
        thread_list = []
        for t in range(min(thread,len(re_filter_dict))): #减少进程开销
            p = Process(target=do_re_filter_loop,
                        args=(
                        parameter_information, re_filtered_reads_whole_bp_dict, ref_length_dict, re_filter_dict, t, min(thread,len(re_filter_dict))))
            thread_list.append(p)
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()
        t3 = time.time()
        print("Time used: {}s".format(round((t3 - t2), 4)))
        print("Whole time: {}s".format(round((t3 - t1), 4)))
    else:
        t3 = time.time()
        print("Time used: {}s".format(round((t3 - t2), 4)))
        print("Whole time: {}s".format(round((t3 - t1), 4)))


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="filter",usage="%(prog)s <-1 -2| -s> <-r> <-o> [options]")
    pars.add_argument("-1", "--data1", dest="data1",
                      help="One end of the paired-end reads, support fastq/fastq.gz/fastq.bz2", metavar="")
    pars.add_argument("-2", "--data2", dest="data2",
                      help="Another end of the paired-end reads, support fastq/fastq.gz/fastq.bz2",
                      metavar="")

    pars.add_argument("-s", "--single", dest="single",
                      help="Single-read, support fastq/fastq.gz/fastq.bz2", metavar="")
    pars.add_argument("-o", "--out", dest="out", help="Specify the result folder",
                      metavar="", required=True)

    pars.add_argument("-k1", "--kmer2", dest="wordsize", help="length of a word-size [default=17]",
                      default=17, type=int, metavar="")
    pars.add_argument("-step_length", metavar="", dest="step", type=int,
                      help="the length of the sliding window on the reads, [default=4]", default=4)  # step length
    pars.add_argument("-t", "--thread", metavar="", dest="thread", help="Thread", type=int,default=1)
    pars.add_argument("-r", "--reference", metavar="", dest="reference", type=str, help="references", required=True)

    pars.add_argument("-d","--data", metavar="", dest="data_size", help="data size", default='all')
    args = pars.parse_args()

    my_filter_main(args)






















