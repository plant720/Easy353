#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/4/14 17:35
# @Author  : xiepulin
# @File    : my_assemble_new.py
# @Software: PyCharm
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/3/31 9:32
# @Author  : xiepulin
# @File    : my_assemble.py
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
from concurrent.futures import ProcessPoolExecutor


'''
获取绝对路径
'''
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

'''
获得目标文件夹下的fasta文件，仅第一层
'''
def get_fasta_file(path):
    path=get_absolute(path)
    files=os.listdir(path)
    suffix_list=["fa","fas","fasta"]
    files_list=[]
    for i in files:
        suffix=i.split(".")[-1].lower()
        file_path= os.path.join(path,i)
        if os.path.isfile(file_path) and suffix in suffix_list :
            files_list.append(file_path)
    return  files_list


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



# 反向互补
def reverse_complement_all(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
# 简化反向互补
def reverse_complement_limit(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

# kmers生成器
def make_kmers(seq, k):
    for i in range(len(seq)-k+1): yield seq[i:i+k]

# 补齐生成器
def forward(seq):
    for x in 'ACGT': yield seq[1:]+x
def reverse(seq):
    for x in 'TGCA': yield x + seq[:-1]

# 创建哈西字典,将reference存入字典

'''
data_structure {"kmer1":[kmer1_count,pos1,gene_name1,gene_name2,gene_name3],
                "kmer2":[kmer2_count,pos2,gene_name1,gene_name2,gene_name3]
                }

{'CGTTTTGACTGTATCGC': [3, 0.0017985611510791366, 'matK'],
'GCGATACAGTCAAAACG': [3, -2.9712230215827335, 'matK']}
'''
def gethashdict(reference, merSize, get_rc = False,pos = False,print_info=False):
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
            refseq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))  #将一个文件下所有seq连在一起
            for j in range(0, len(refseq)-merSize+1):
                temp_list, kmer = [], refseq[j:j+merSize]
                if kmer in kmer_dict:
                    temp_list = kmer_dict[kmer]
                    temp_list[0] += 1
                    if pos:
                        temp_list[1] += (j+1)/len(refseq)  #(j+1)/len(refseq) 记录位置信息，如果反复出现，我们可以理解为在求一个位置的均值 （估计值）
                    if file_name not in temp_list:
                        temp_list.append(file_name)
                else:
                    temp_list = [1, (j+1)/len(refseq), file_name] if pos else [1, 0, file_name] # 记录位置信息
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
                        if pos:
                            temp_list[1] -= (j+1)/len(refseq)
                        if file_name not in temp_list: temp_list.append(file_name)
                    else:
                        temp_list = [1, -(j+1)/len(refseq), file_name] if pos else [1, 0, file_name]
                        kmer_dict[sys.intern(kmer)] = temp_list  #sys,inter利用了地址信息，保证相同元素出现在同一地址
            if print_info:
                print('Memory:', round(sys.getsizeof(kmer_dict)/1024/1024, 4), 'MB, Number of Genes:', gene_number, end='\r',flush=True)
        infile.close()
    return  kmer_dict

'''
assemble的hashdict  记录kmercount频次，记录kmer在ref的平均位置信息
'''
def make_assemble_hashdict(filtered_out, merSize,ref_dict,get_rc = False):
    files=get_file_list(filtered_out)
    if files==[]:
        return 0
    kmer_dict=defaultdict(list)
    kmer_count=0
    for file in files:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            refseq_r = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))  #过滤器，返回迭代器对象
            if get_rc:
                for refseq in [refseq_r, reverse_complement_limit(refseq_r)]:
                    kmer_count += len(refseq) - merSize + 1
                    for j in range(0, len(refseq) - merSize + 1):
                        temp_list, temp_len, kmer = [], [0], refseq[j:j + merSize]
                        if kmer in kmer_dict:
                            temp_list = kmer_dict[kmer]
                            temp_list[0] += 1
                        else:
                            temp_list = [1, temp_len]
                            kmer_dict[kmer] = temp_list
                        if kmer in ref_dict:
                            temp_len[0] = 1 + ref_dict[kmer][1] / ref_dict[kmer][0] if ref_dict[kmer][1] < 0 else ref_dict[kmer][1] / ref_dict[kmer][0]

                        # ref_dict[kmer][1] / ref_dict[kmer][0]  累计位置信息/kmercount = 平均位置信息  之所以用1+，是因为考虑了反向互补的情况

            else:
                for refseq in [refseq_r]:
                    kmer_count += len(refseq) - merSize + 1
                    for j in range(0, len(refseq)-merSize+1):
                        temp_list, temp_len,  kmer = [], [0], refseq[j:j+merSize]
                        if kmer in kmer_dict:
                            temp_list = kmer_dict[kmer]
                            temp_list[0] += 1
                        else:
                            temp_list = [1, temp_len]
                            kmer_dict[kmer] = temp_list
                        if kmer in ref_dict:
                            temp_len[0] = 1 + ref_dict[kmer][1] / ref_dict[kmer][0] if ref_dict[kmer][1] < 0 else ref_dict[kmer][1] / ref_dict[kmer][0]

        infile.close()

    return  kmer_dict,kmer_count



def get_contig_forward(_dict, seed, iteration = 1000 ,weight = 0.5):
    '''
    :param _dict:  filtered_out哈希字典
    :param seed:   seed
    :param iteration:  迭代次数 之前指最大递归次数
    :return:
    '''

    '''
    temp_list:  存储用于组装的kmer,每一次弹出最后的一个kmer,用于延申
    reads_list: 存储用于组装的kmer
    stack_list: 栈 将Kmer记录为节点，遍历树的方法，每一次到叶子节点就回溯到上一个分歧节点，直到栈为空（出现分歧的时候记录）
    best_kmc:  记录最大累计Kmercount
    cur_kmc: 记录当前累计kmercount
    best_pos 记录最佳位置 ，用于contig组装scaffold
    cur_pos   记录当前位置
    best_seq:记录最佳路径
    cur_seq:记录当前路径
    
    核心，kmer的组装是按照权重组装的，权重由reads的kmercout 和reads在ref上的位置决定
    
    _dict[i][0] ** get_dis(_pos, _dict[i][1][0], weight)  没有位置信息，将开根号处理（默认0.5）;    核心思想距离越近，kmercount越大的权重越大
    
    '''
    temp_list, reads_list, stack_list = [seed], [seed], []
    best_kmc, cur_kmc, best_pos, cur_pos, best_seq, cur_seq = [], [], [], [], [], []
    _pos, node_distance = 0, 0
    while True:
        node = sorted([(i, _dict[i][1][0], _dict[i][0] ** get_dis(_pos, _dict[i][1][0], weight)) for i in forward(temp_list[-1]) if i in _dict], key=lambda _: _[2], reverse=True)
        # print(node)
        #node  [ "kmer",cur_pos, weight]  [('GTTTTGACTGTATCGCA', 0.001199040767386091, 24.995845787996004)]
        while node:
            if node[0][0] in temp_list:
                node.pop(0)
            else:
                break
        if not node:
            iteration -= 1
            if sum(cur_kmc) > sum(best_kmc):
                best_kmc, best_seq, best_pos = cur_kmc.copy(), cur_seq.copy(), cur_pos.copy()
            [(temp_list.pop(), cur_kmc.pop(), cur_seq.pop(), cur_pos.pop()) for _ in range(node_distance)]
            if stack_list:
                node, node_distance, _pos = stack_list.pop()
            else:
                break
        if len(node) >= 2:
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0:
            _pos = node[0][1]
        temp_list.append(node[0][0])
        reads_list.append(node[0][0])
        cur_pos.append(node[0][1])             #要么有位置信息，要么未出现在ref上则为0  有位置信息的会被优先组装，随后在原有基础上向两边延申
        cur_kmc.append(node[0][2])
        cur_seq.append(node[0][0][-1])
        node_distance += 1
        if not iteration: break
    return best_seq, reads_list, best_kmc, best_pos


def get_contig(_dict, seed, iteration = 1000 ,weight = 0.5):
    right, reads_list_1, k_1, pos_1 = get_contig_forward(_dict, seed, iteration, weight) #best_seq, reads_list, best_kmc, best_pos
    left, reads_list_2, k_2, pos_2 = get_contig_forward(_dict, reverse_complement_limit(seed), iteration, weight)

    # print("pos_1:{}".format(pos_1))   #[0.1,0.2,0.3...0.9, 0 , 0, 0, 0 .0 ,0 ,0 ]
    # print("pos_2:{}".format(pos_2))   #[0 , 0, 0, 0 .0 ,0 ,0 ]

    _pos = [x for x in pos_1 + pos_2 if x > 0]
    min_pos, max_pos = 0, 1
    if _pos: min_pos, max_pos = min(_pos), max(_pos)
    return reverse_complement_limit(''.join(left)) + seed + ''.join(right), min_pos, max_pos, reads_list_1 + reads_list_2


def get_scaffold(_dict, seed_list, limit, ref_length, iteration=1000, weight=0.5):
    mid_seed, min_dis = seed_list[0], 1
    # for seed in seed_list:
    #     if abs(seed[2] - 0.5) < min_dis:
    #         mid_seed, min_dis = seed, abs(seed[2] - 0.5)

    contig, min_pos, max_pos, read_list = get_contig(_dict, mid_seed[0], iteration, weight)

    for i in read_list:
        if i in _dict[i]: del _dict[i]
    contigs = [(contig, min_pos, max_pos)]
    # print(seed_list) #[('TTTGGTGGAGTTGTTGGTG', 92, 79.03477756740409)] ["kmer",kmercount.pos_sum累计位置]



    #左侧拼接 left
    new_seed = []
    for x in seed_list:
        if x[2] / x[1] < min_pos and x[1] > limit:
            new_seed = x
    while new_seed:
        #将使用过的reads剔除
        contig, _min_pos, _max_pos, read_list = get_contig(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]:
                del _dict[i]

        # print(new_seed[2] / new_seed[1],"1")
        # print(min_pos,max_pos,"pos")
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos:
            break
        # print(_min_pos, _max_pos,contig,"left")
        if _max_pos < contigs[0][2]:   #新contig完全在左侧  [(contig, min_pos, max_pos)]
            contigs.insert(0, (contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] < _min_pos and x[1] > limit:
                new_seed = x

    #右侧拼接
    new_seed = []
    for x in seed_list:
        if x[2] / x[1] > max_pos and x[1] > limit:
            new_seed = x
    while new_seed:
        contig, _min_pos, _max_pos, read_list = get_contig(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]:
                del _dict[i]
        # print(new_seed[2] / new_seed[1],"1")
        # print(min_pos,max_pos,"pos")
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos:
            break
        # print(_min_pos, _max_pos, contig,"right")
        if _min_pos > contigs[-1][1]:   # 完全在右侧 [(contig, min_pos, max_pos)]
            contigs.append((contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] > _max_pos and x[1] > limit:
                new_seed = x

    scaffold = contigs[0][0]
    for x in range(1, len(contigs)):
        scaffold += max(2, int((contigs[x][1] - contigs[x - 1][2]) * ref_length)) * 'N' + contigs[x][0]  #contig 从左往右组装， 下一个最小contig最小位置减去上一个contig最大位置 能估算出gap空洞长度
    return scaffold









'''
利用交集确定种子，利用动态limit过滤低质量reads
'''
def get_seed(ref_dict,filtered_dict,ref_length,limit_count):
    #利用交集确定种子
    for i in ref_dict:
        if i not in filtered_dict:
            ref_dict[i][0] = 0
    #动态limit
    if limit_count < 0:
        limit_count = dynamic_limit_v2(filtered_dict, ref_length, 2)
    #剔除低质量reads: 频次太低且未在ref中出现，注意频次很高未在ref中的reads是保留了的
    if limit_count > 0:
        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []
    seed_list = [(x, ref_dict[x][0]) for x in ref_dict if ref_dict[x][0] > 0]
    list.sort(seed_list, key=lambda x: x[1], reverse=True)
    return  seed_list


def get_dis(_pos, new_pos, weight = 0.5):
    return 1 - abs(_pos - new_pos) if (_pos and new_pos) else weight


def dynamic_limit_v1(_dict, smoother_level=4, smoother_size=64):
    # 计算kmercount 频次分布图,默认记录kmer最大频次为256
    count_list = [0] * smoother_size
    for x in _dict:
        if _dict[x][0] <= smoother_size:
            count_list[_dict[x][0] - 1] += 1

    # print(count_list)
    # 平滑器
    for i in range(smoother_level):
        for x in range(1, smoother_size - smoother_level + i):
            if count_list[x + smoother_level - 1 - i] < count_list[x - 1] and count_list[x] < count_list[x + smoother_level - i]:

                return x + 1
    return 2


#F0=F1+F2+length*rate rate取值为2的时候，比较宽松
def dynamic_limit_v2(_dict, ref_length, N, ref_rate=1.5, list_size=256):
    '''
    :param _dict:   filtered dict
    :param ref_length:   参考序列
    :param N:       filtered_kmer_count
    :param ref_rate:   参考序列倍率
    :param list_size:  记录kmer频次分布
    :return: (limit,deepth)
    理论支撑：
    F0-f1-f2= rate * ref_length*2    2代表2倍，因为计算了反向互补
    (kmercount_all-低频kmercount_sum)  / ref_length  ~= deepth
    '''
     #计算kmercount 频次分布图,默认记录kmer最大频次为256
    count_list = [0] * list_size
    for x in _dict:
        if _dict[x][0] <= list_size:
            count_list[_dict[x][0] - 1] += 1   #[46516, 5624, 1334, 376, 630, 108, 72, 44, 98, 94, 2,2,2,2,2,2,0,0,0,0,2]



    # 计算器
    F0, sum_f, sum_k = len(_dict), 0, 0
    for x in range(list_size):
        sum_f += count_list[x]
        sum_k += count_list[x] * (x + 1)
        if (F0 - sum_f) / 2 < ref_length * ref_rate:
            # print(N, sum_k, ref_length, ref_rate)   #2306410 63270 1514 2
            return max(2, x),  round((N - sum_k) / ref_length / ref_rate / 2, 2) #返回Limit,deepth

    return 2, round((N - sum_k) / ref_length / ref_rate / 2, 2)  #返回Limit,deepth


def dynamic_weight(_dict):  #匹配上的reads越多， 值越小。
    temp_count = 0
    for i in _dict:
        if _dict[i][1][0]>0:  #_dict[i][1][0] 位置信息 大于0说明该条reads匹配到了ref
            temp_count += 1
    return 1 - temp_count/len(_dict) #temp_count/len(_dict)代表相似度  1-相似度=差异度  weight=kmercount**weight 所以，差异度越大，ref不可信，则需要越依赖kmercount的大小


def do_assemble(ref_path,filtered_path,assemble_kmer,limit_count,limit_min_length,limit_max_length,ref_length_dict,assembled_path,change_times=32,scaffold_or_not=True):


    t1=time.time()
    gene_name=get_basename(ref_path)
    ref_dict = gethashdict(ref_path, assemble_kmer, get_rc=True, pos=True, print_info=False)
    filtered_dict, filtered_kmer_count = make_assemble_hashdict(filtered_path, assemble_kmer, ref_dict, get_rc=True)
    ref_length=ref_length_dict[gene_name]
    print(filtered_kmer_count, "new_kmer_sum")


    # 利用交集确定种子
    for i in ref_dict:
        if i not in filtered_dict:
            ref_dict[i][0] = 0
    # 动态limit
    if limit_count =="auto":
        limit_count_v2, deepth = dynamic_limit_v2(filtered_dict, ref_length, filtered_kmer_count, 2)
        print("limit_count_v2:{}".format(limit_count_v2))
        if deepth > 20:  # 对于测序深度较深的序列，可以使用较为严格的limit
            limit_count_v1 = dynamic_limit_v1(filtered_dict, 8, 64)
            limit_count = min(limit_count_v2,limit_count_v1)
            print("limit_count:{0},limit_count_v1:{1},limit_count_v2:{2},deepth:{3}".format(limit_count,limit_count_v1,limit_count_v2,deepth))

    # 剔除低质量reads: 频次太低且未在ref中出现，注意频次很高未在ref中的reads是保留了的
    if limit_count > 0:
        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []

    print(len(filtered_dict), "new_kmer_sum")
    # seed_list = [(x, ref_dict[x][0]) for x in ref_dict if ref_dict[x][0] > 0]
    seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if ref_dict[x][0] > 0 and ref_dict[x][1] > 0] #新增位置信息
    if seed_list == []:
        return 0
    list.sort(seed_list, key=lambda x: x[1], reverse=True)





    # 利用相似度计算权重 filtered reads 与ref 的大致距离
    cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5  # 0.52
    print(cur_weight, "cur_weight")

    '''
    best_seed,best_contig,best_length 设定第一次默认最佳
    '''
    best_seed=seed_list[0]
    best_contig = get_contig(filtered_dict, best_seed[0], iteration=1000, weight=cur_weight)[0]  # 速度与允许迭代次数高度正相关，python最大递归1000次，我们改写为循环后，可突破上限
    best_length=len(best_contig)

    ref_length=ref_length_dict[gene_name]
    contig_path = os.path.join(assembled_path, "contig", gene_name + ".fasta")
    short_contig_path=os.path.join(assembled_path, "short_contig", gene_name + ".fasta")
    scaffold_path = os.path.join(assembled_path, "scaffold", gene_name + ".fasta")


    ref_length_min=ref_length*limit_min_length
    ref_length_max=ref_length*limit_max_length


    write_contig=False
    if best_length < ref_length_min:
        #更换seed
        change_time=0
        last_seed = seed_list[0]
        for seed in seed_list[1:]:
            if last_seed[0][1:] == seed[0][:-1] or last_seed[0][2:] == seed[0][:-2]: #移动1~2bp
                last_seed = seed
                continue
            change_time+=1
            last_seed = seed
            new_contig=get_contig(filtered_dict, best_seed[0], iteration=1000, weight=cur_weight)[0]
            new_length=len(new_contig)

            # print("times:{0} length:{1} contig:{2}".format(change_time,new_length,new_contig)) #for test
            print("{0}: use another seed: {1} ... {2:<2} ".format(gene_name,last_seed[0][0:8],change_time))
            if new_length>=ref_length_min:
                best_seed=last_seed
                best_contig=new_contig
                best_length=new_length
                write_contig=True
                break
            elif change_time>=change_times: #超过最大更换种子数目
                #种子为未达标情况下，最佳seed
                break
            else: #未达标情况下的最好结果
                if new_length>best_length:
                    best_contig=new_contig
                    best_length=new_length
                    best_seed=last_seed


    #对于过长的处理，不更换seed，改为增加limit
    elif best_length>ref_length_max:
        if limit_count >=32:
            limit_count_max=limit_count+32
        else:
            limit_count_max=33
        # 加大limit
        for temp_limit in range(limit_count+1,limit_count_max):
            # 对filtered_dict再次精简
            _filter = [x for x in filtered_dict if filtered_dict[x][0] <= temp_limit]
            for x in _filter:
                if filtered_dict[x][1][0] == 0:
                    del filtered_dict[x]
            _filter = []
            # 此处更新权重影响很小，不更新
            # cur_weight = assemble.dynamic_weight(filtered_dict)
            contig=get_contig(filtered_dict, best_seed[0], iteration=1000, weight=cur_weight)[0]
            best_contig = contig
            best_length=len(best_contig)
            best_seed=best_seed

            # print("{0}: add limit ".format(gene_name))
            if best_length < ref_length_max:
                break
        write_contig = True
    else:
        write_contig = True



    if write_contig:
        print("{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8], best_length,
                                                                              ref_length))
        with open(contig_path, "w") as f:
                f.write('>'+gene_name+'_contig_k'+str(assemble_kmer)+'_'+str(best_length)+'\n')
                f.write(best_contig+'\n')
    else:
        print("{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8], best_length,
                                                                              ref_length))
        with open(short_contig_path, "w") as f:
                f.write('>' +gene_name+ '_short_contig_k' + str(assemble_kmer) + '_' + str(best_length) + '\n')
                f.write(best_contig + '\n')



    t2=time.time()
    if scaffold_or_not:
        seed_list.insert(0, best_seed)
        scaffold = get_scaffold(filtered_dict, seed_list, limit_count, ref_length, weight=cur_weight)
        with open(scaffold_path, 'w') as out:
            scaffold_length = len(scaffold.replace("N", ""))
            out.write('>'+gene_name+'_scaffold_k' +str(assemble_kmer)+'_'+ str(scaffold_length) + '\n')
            out.write(scaffold + '\n')
            if not write_contig:
                if scaffold_length > limit_min_length:
                    print("{}: try scaffold. scaffold length: {}".format(gene_name, scaffold_length))
                else:
                    pass

    t3=time.time()
    print(t3-t1,t2-t1,t3-t2)

    if is_exist(contig_path) or is_exist(scaffold_path):
        if is_exist(short_contig_path):
            os.remove(short_contig_path)




def do_assemble_parallel(ref_all_path,filtered_all_path,assemble_kmer,limit_count,limit_min_length, limit_max_length,ref_length_dict,assembled_path,thread,change_times=32,scaffold_or_not = True):



    files=get_fasta_file(filtered_all_path)
    task_pool=[]
    results=[]

    if not os.path.isdir(assembled_path):
        os.mkdir(assembled_path)
    short_contig_path=os.path.join(assembled_path, "short_contig")
    contig_path=os.path.join(assembled_path,"contig")
    scaffold_path=os.path.join(assembled_path,"scaffold")

    if not os.path.isdir(short_contig_path):
        os.mkdir(short_contig_path)
    if not os.path.isdir(contig_path):
        os.mkdir(contig_path)
    if not os.path.isdir(scaffold_path) and scaffold_or_not:
        os.mkdir(scaffold_path)

    executor=ProcessPoolExecutor(max_workers=thread)
    for i in files:
        name=get_basename(i)
        ref_path=os.path.join(ref_all_path,name+".fasta")
        task_pool.append(executor.submit(do_assemble,ref_path,i,assemble_kmer,limit_count,limit_min_length, limit_max_length,ref_length_dict,assembled_path,change_times,scaffold_or_not))

    for i in task_pool:
        results.append(i.result())
    executor.shutdown()




def my_assemble_main(args):
    # 初始化操作
    thread = args.thread
    kmer = args.kmer
    assembled_path= args.out
    reference = args.reference
    filtered_path=args.input
    limit_count=args.limit_count
    limit_min_length=args.limit_min_length
    limit_max_length = args.limit_max_length
    change_times=args.change_times
    scaffold_or_not=args.scaffold


    True_set=["yes","y","Yes","true","True","1",1]
    if scaffold_or_not in True_set:
        scaffold_or_not =True
    else:
        scaffold_or_not =False

    print('======================== Assemble =========================')
    t1 = time.time()
    ref_length_dict = {}
    ref_length_dict=get_ref_info(reference,ref_length_dict)
    do_assemble_parallel(reference, filtered_path,kmer, limit_count, limit_min_length, limit_max_length,
                ref_length_dict, assembled_path,thread,change_times, scaffold_or_not)

    t2 = time.time()
    print("time used: {}".format(t2-t1))


# if __name__ == '__main__':
#     t1=time.time()
#     # ref_path = r"example\matk_ref"
#     # filtered_path = r"example\matk_filtered"
#     # assembled_path=r"example\matk_assembled"
#     #
#     ref_path = r"example\5264_ref\5264.fasta"
#     filtered_path = r"example\5264_filtered\5264.fasta"
#     assembled_path = r"example\5264_assembled"
#
#
#     assemble_kmer = 17
#     limit= -1
#     limit_length = 1
#     thread=4
#     change_times = 32
#     # value = 1514
#     value=1668
#     ref_length_dict={}
#     ref_length_dict=get_ref_info(ref_path,ref_length_dict)
#     #分开写
#     ref_dict = gethashdict(ref_path, assemble_kmer, get_rc=True, pos=True, print_info=False)
#     filtered_dict,filtered_kmer_count = make_assemble_hashdict(filtered_path, assemble_kmer, ref_dict, get_rc=True)
#     # print(len(filtered_dict))
#     print(filtered_kmer_count)
#
#
#     for i in ref_dict:
#         if i not in filtered_dict:
#             ref_dict[i][0] = 0
#     if limit < 0:
#         limit,deepth = dynamic_limit_v2(filtered_dict, value,filtered_kmer_count, 2)
#         if deepth > 20: #对于测序深度较深的序列，可以使用较为严格的limit
#             limit = max(limit, dynamic_limit_v1(filtered_dict, 8, 64))
#         print('\n', limit,deepth)
#
#     limit =3
#     if limit > 0:
#         _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit]
#         for x in _filter:
#             if filtered_dict[x][1][0] == 0:  #剔除无效reads  低于limit且未在ref中出现的reads
#                 del filtered_dict[x]
#         _filter = []
#     # print(len(filtered_dict)) #59840 -- 6446
#     # print(filtered_dict)
#     cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5    #0.52
#     # print(cur_weight)
#     seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if ref_dict[x][0] > 0 and ref_dict[x][1]>0 ]
#     list.sort(seed_list, key=lambda x: x[1], reverse=True)  #seed_list ef与reads取交集之后 按照kmercount排序  [("kmer1",kmercount,pos),("kmer1",kmercount,pos)]
#     t1 = time.time()
#
#     last_seed = seed_list[0][0]
#     contig = get_contig(filtered_dict, last_seed, weight=cur_weight)[0]
#     scaffold = get_scaffold(filtered_dict, seed_list, limit, value)
#     print("\n>contig_" + str(len(contig)))
#     print(contig)
#     print("\n>scaffold_" + str(len(contig)))
#     print(scaffold)
#     t2 = time.time()
#     print("time:",t2-t1)


# if __name__ == '__main__':
#     t1=time.time()

    # # ref_path = r"example\matk_ref\matK.fasta"
    # # filtered_path = r"example\matk_filtered\matK.fasta"
    # # assembled_path = r"example\matk_assembled"
    #
    # # ref_path = r"example\trnk_ref\trnK-UUU.fasta"
    # # filtered_path = r"example\trnk_filtered\trnK-UUU.fasta"
    # # assembled_path = r"example\trnk_assembled"
    #
    # ref_path = r"example\5264_ref\5264.fasta"
    # filtered_path = r"example\5264_filtered\5264.fasta"
    # assembled_path = r"example\5264_assembled"
    #
    # ref_path_all=r"example\cp_gene"
    # filtered_path_all=r"example\filtered_out"
    # assembled_path_all=r"example\cp_assembled"
    #
    # limit_count = -1
    # assemble_kmer = 19
    # limit_min_length=0.5
    # limit_max_length=2
    # thread=4
    # change_times = 32
    # scaffold_or_not=True
    #
    # ref_length_dict={}
    # ref_length_dict=get_ref_info(ref_path,ref_length_dict)
    #
    # ref_length_dict_all={}
    # ref_length_dict_all=get_ref_info(ref_path_all,ref_length_dict_all)
    #
    # do_assemble(ref_path,filtered_path,assemble_kmer,limit_count,limit_min_length,limit_max_length,ref_length_dict,assembled_path,change_times=32,scaffold_or_not=True)
    # # do_assemble_parallel(ref_path_all, filtered_path_all, assemble_kmer, limit_count,limit_min_length,limit_max_length, ref_length_dict_all, assembled_path_all,thread,  change_times,scaffold_or_not)
    #



if __name__ == '__main__':
    #-i -r 中的文件名必须是xx.fasta 同名
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="assemble",
                                   usage="%(prog)s <-i> <-r> <-o> [options]")

    pars.add_argument("-o", "--out", dest="out", help="Specify the result folder <dir>",
                      metavar="", required=True)

    pars.add_argument("-k2", "--kmer2", dest="kmer", help="length of kmer[default=31]",
                      default=31, type=int, metavar="")
    pars.add_argument("-t", "--thread", metavar="", dest="thread", help="Thread", type=int, default=1)
    pars.add_argument("-r", "--reference", metavar="", dest="reference", type=str, help="references  <dir>", required=True)
    pars.add_argument("-i", "--input", metavar="", dest="input", type=str, help="filtered reads  <dir>", required=True)

    pars.add_argument('-limit_count', metavar='',dest='limit_count', help='''limit of kmer count [default=2]''', required=False,
                      default=2)
    pars.add_argument('-limit_min_length', metavar='', dest='limit_min_length',type=float, help='''limit of contig length''',
                      required=False, default=0.5)
    pars.add_argument('-limit_max_length', metavar='',dest='limit_max_length', type=float, help='''limit of contig length''',
                      required=False, default=2)
    pars.add_argument("-change_seed", metavar="", dest="change_times",type=int, help='''times of changing seed [default=32]''', required=False,
                      default=32)
    pars.add_argument('-scaffold', metavar="",dest="scaffold",type=str,help='''make scaffold''', default=True)
    args = pars.parse_args()

    my_assemble_main(args)









