#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 4:59 下午
# @Author     : zzhen
# @File       : assembly.py
# @Software   : PyCharm
# @Description:
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import os
from shutil import rmtree
import argparse
import collections

# import my python file
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from Easy353Lib.utils import make_ref_kmer_dict, get_seq_avg_length, reverse_complement_limit, \
    get_reads_info, get_file_list, log


# 获取list的中位数
def get_median(lst: list) -> float:
    lst.sort()
    half = len(lst) // 2
    return (lst[half] + lst[~half]) / 2


# kmers生成器
def make_kmers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]


# 补齐生成器
def forward(seq):
    for x in 'ACGT':
        yield seq[1:] + x


def reverse(seq):
    for x in 'TGCA':
        yield x + seq[:-1]


# 将测序数据中的read进行处理 并返回处理后的read_kmer_dict
# read_kmer_dict用于记录过滤后的read拆分后的kmer key:read_kmer value:[read_count,read_pos]
# _read_file_list_:用于拼接的测序文件列表 _kmer_size_ _ref_kmer_dict_:参考序列的kmer字典
def make_assemble_hash_dict(_read_file_: str, _ref_file_: str, _kmer_size_: int,
                            _ref_reverse_complement_: bool = False, _read_reverse_complement_: bool = False,
                            _pos_: bool = True, _print_: bool = True) -> tuple:
    # 从参考文件中获取种子字典
    _ref_kmer_dict_ = make_ref_kmer_dict(_ref_file_, _kmer_size_, _ref_reverse_complement_, _pos_,
                                         _print_=False)
    # 获取参考序列平均长度
    _ref_avg_length = get_seq_avg_length(_ref_file_)
    read_kmer_dict = collections.defaultdict(list)
    infile = open(_read_file_, 'r', encoding='utf-8', errors='ignore')
    infile.readline()
    for line in infile:
        temp_str = []
        while line and line[0] != '>':
            temp_str.append(line)
            line = infile.readline()
        _read_seq_ = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
        # 该处由于已经是成对的reads，所以不应该对reads进行返回反向互补的操作
        reads_seqs = [_read_seq_]
        if _read_reverse_complement_:
            reads_seqs.append(reverse_complement_limit(_read_seq_))
        for _seq_ in reads_seqs:
            # 用于标识当前read产生的kmer是否在参考序列中，如果一个kmer在参考序列中，则flag置为true,记录位置信息，
            # 并将后续kmer也添加在参考序列中
            kmer_flag = False
            # 分别代表距离上一个被记录的kmer的距离pos_int，当前kmer的距离temp_pos
            pos_list = [1, 0]
            for j in range(0, len(_seq_) - _kmer_size_ + 1):
                kmer_info_list, read_kmer = [], _seq_[j:j + _kmer_size_]
                if read_kmer in read_kmer_dict:
                    read_kmer_dict[read_kmer][0] += 1
                else:
                    # 只要一条read的kmer出现在kmer_dict中，就统一记录kmer的信息
                    # 在ref_dict中的位置总信息/ref_dict中出现的总次数
                    if read_kmer in _ref_kmer_dict_:
                        # 在参考序列的方向互补序列上
                        if _ref_kmer_dict_[read_kmer][1] < 0:
                            pos_list[1] = 1 + _ref_kmer_dict_[read_kmer][1] / _ref_kmer_dict_[read_kmer][0]
                        else:
                            pos_list[1] = _ref_kmer_dict_[read_kmer][1] / _ref_kmer_dict_[read_kmer][0]
                        kmer_flag = True
                    else:
                        if kmer_flag:
                            pos_list[1] += pos_list[0] / _ref_avg_length
                    # 当一个kmer既不在参考序列中，kmer_flag也不为true时，则位置信息为0
                    kmer_info_list = [1, pos_list[1]]
                    read_kmer_dict[read_kmer] = kmer_info_list
    infile.close()
    return read_kmer_dict, _ref_kmer_dict_


# temp_list:  存储用于组装的kmer,每一次弹出最后的一个kmer,用于延申
# reads_list: 存储用于组装的kmer
# stack_list: 栈 将Kmer记录为节点，遍历树的方法，每一次到叶子节点就回溯到上一个分歧节点，直到栈为空（出现分歧的时候记录）
# best_kmc: 记录最大累计节点权重值
# cur_kmc: 记录当前累计节点权重值
# best_pos: 记录最佳位置 ，用于contig组装scaffold
# cur_pos: 记录当前位置
# best_seq: 记录最佳路径
# cur_seq: 记录当前路径
def assemble_contig_forward(_read_kmer_dict_: dict, _seed_: str, iteration: int = 1000) -> tuple:
    temp_list, stack_list = [_seed_], []
    best_kmer_weight, cur_kmer_weight, best_seq, cur_seq, best_pos, cur_pos = [], [], [], [], [], []
    # 初始化一个set 每一个元素是一个元组 第一个元素是拼接过程中使用过的kmer，第二个元素是kmer的数量
    used_kmer_info = {(_seed_, _read_kmer_dict_.get(_seed_)[0])}
    _pos, node_distance = 0, 0
    while True:
        # 动态生成下次迭代的kmer 并根据权重进行排列
        node = sorted(
            [(i, _read_kmer_dict_[i][1], count_pos(_read_kmer_dict_[i][0], _pos, _read_kmer_dict_[i][1]),
              _read_kmer_dict_[i][0]) for i
             in forward(temp_list[-1])
             if i in _read_kmer_dict_], key=lambda _: _[2], reverse=True)
        while node:
            if node[0][0] in temp_list:
                node.pop(0)
            else:
                break
        if not node:
            iteration -= 1
            if sum(cur_kmer_weight) > sum(best_kmer_weight):
                best_kmer_weight, best_seq, best_pos = cur_kmer_weight.copy(), cur_seq.copy(), cur_pos.copy()
            [(temp_list.pop(), cur_kmer_weight.pop(), cur_seq.pop(), cur_pos.pop()) for _ in range(node_distance)]
            if stack_list:
                node, node_distance, _pos = stack_list.pop()
            else:
                break
        if len(node) >= 2:
            # 放入stack_list，可用于下次迭代时进行回溯
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0:
            _pos = node[0][1]
        temp_list.append(node[0][0])
        cur_pos.append(node[0][1])
        cur_kmer_weight.append(node[0][2])
        cur_seq.append(node[0][0][-1])
        # 在计算使用过的kmer数量时，回溯的kmer数量不去除
        used_kmer_info.add((node[0][0], node[0][3]))
        node_distance += 1
        if not iteration:
            break
    # 返回拼接出的contig,最终的位置信息, 使用的reads信息 和 使用的kmer数量
    reads_list = [i[0] for i in used_kmer_info]
    cur_used_kmer_count = [i[1] for i in used_kmer_info]
    return best_seq, best_pos, reads_list, sum(cur_used_kmer_count)


# 用于将read的出现次数和位置信息进行整合
def count_pos(count, cur_pos, average_pos):
    if cur_pos and average_pos:
        return count ** float(1 - abs(cur_pos - average_pos))
    return count ** 0.5


def assemble_read_for_contig(_read_kmer_dict_: dict, _seed_: str, iteration: int = 1000):
    right, pos_1, reads_list_1, kmer_count_1 = \
        assemble_contig_forward(_read_kmer_dict_, _seed_, iteration)
    left, pos_2, reads_list_2, kmer_count_2 = \
        assemble_contig_forward(_read_kmer_dict_, reverse_complement_limit(_seed_), iteration)
    # 返回一个元组 (拼接的contig,reads_list包含所有被访问过的reads)
    _pos_all_ = [x for x in pos_1 + pos_2 if x > 0]
    min_pos, max_pos = 0, 1
    if _pos_all_:
        min_pos, max_pos = min(_pos_all_), max(_pos_all_)
    return (reverse_complement_limit(''.join(left)) + _seed_ + ''.join(right),
            min_pos, max_pos, reads_list_1 + reads_list_2, len(set(reads_list_1 + reads_list_2)),
            kmer_count_1 + kmer_count_2)


def get_scaffold(_read_kmer_dict_: dict, _seed_: list, _kmer_limit_count_: int, gene_avg_len: int,
                 iteration: int = 1000) -> str:
    # 先拼接产生一个contig
    mid_seed, min_dis = _seed_[0], 1
    # 返回一个元组 (拼接的contig,reads_list包含所有被访问过的reads)
    contig, min_pos, max_pos, read_list, _, _ = assemble_read_for_contig(_read_kmer_dict_, mid_seed[0], iteration)
    for i in read_list:
        if i in _read_kmer_dict_:
            del _read_kmer_dict_[i]
    contigs = [(contig, min_pos, max_pos)]
    _left_seed_ = tuple()

    # _seed_中的元素是一个元组 (_kmer_,_count_,avg_pos)
    for x in _seed_[1:]:
        # 选取出现在该contig左边出现的seed用于拼接contig
        if x[2] < min_pos and x[1] > _kmer_limit_count_:
            _left_seed_ = x
    while _left_seed_:
        contig, _min_pos, _max_pos, read_list, _, _ = assemble_read_for_contig(_read_kmer_dict_, _left_seed_[0],
                                                                               iteration)
        for i in read_list:
            if i in _read_kmer_dict_:
                del _read_kmer_dict_[i]
        if _left_seed_[2] < _min_pos or _left_seed_[2] > _max_pos:
            break
        # print(_min_pos, _max_pos,contig)
        if _max_pos < contigs[0][2]:
            contigs.insert(0, (contig, _min_pos, _max_pos))
        else:
            break
        _left_seed_ = tuple()
        for x in _seed_:
            if x[2] < _min_pos and x[1] > _kmer_limit_count_:
                _left_seed_ = x
    _right_seed_ = tuple()
    for x in _seed_:
        # 选取出现在该contig右边的seed用于拼接contig
        if x[2] > max_pos and x[1] > _kmer_limit_count_:
            _right_seed_ = x
    while _right_seed_:
        contig, _min_pos, _max_pos, read_list, _, _ = assemble_read_for_contig(_read_kmer_dict_, _right_seed_[0],
                                                                               iteration)
        for i in read_list:
            if i in _read_kmer_dict_:
                del _read_kmer_dict_[i]
        if _right_seed_[2] < _min_pos or _right_seed_[2] > _max_pos:
            break
        # print(_min_pos, _max_pos, contig)
        if _min_pos > contigs[-1][1]:
            contigs.append((contig, _min_pos, _max_pos))
        else:
            break
        _right_seed_ = tuple()
        for x in _seed_:
            if x[2] > _max_pos and x[1] > _kmer_limit_count_:
                _right_seed_ = x
    scaffold = contigs[0][0]
    for x in range(1, len(contigs)):
        scaffold += max(2, int((contigs[x][1] - contigs[x - 1][2]) * gene_avg_len)) * 'N' + contigs[x][0]
    return scaffold


# _iteration_用于设定迭代次数 在进行assembly时,read是一定需要进行反向互补的 _ref_reverse_complement_
# 在进行assembly时，参考序列需要进行反向互补 这是因为组装contig时选取的种子是反向互补的，当参考序列不进行反向互补时，就会报错
def assemble_one_file(_reads_file_: str, _out_dir_: str, _ref_file_: str,
                      _assemble_kmer_size_: int, _ref_reverse_complement_: bool = True,
                      _pos_: bool = True, _change_seed_: int = 1000, _kmer_limit_count_: int = 2,
                      _min_percent_length_: float = 0.5, _max_percent_length_: float = 1.5,
                      _iteration_: int = 1000, _print_: bool = True) -> dict:
    assemble_gene_info_dict = collections.defaultdict(dict)
    # 设定输出文件名
    gene_name = os.path.basename(_ref_file_).split(".")[0]
    target_file_name = gene_name + ".target.fasta"
    unrecovered_dir = os.path.join(_out_dir_, "unrecovered_genes")

    # 根据_ref_file_获取gene的平均长度和ref中的序列数量
    ref_seq_count, ref_seq_length = get_reads_info(_ref_file_)
    gene_avg_len = round(ref_seq_length / ref_seq_count, 2)

    _contig_path_ = os.path.join(_out_dir_, target_file_name)
    _short_contig_path_ = os.path.join(unrecovered_dir, target_file_name)

    assemble_success_flag = False
    is_normal_contig = False
    gene_info = {"assemble_success_flag": assemble_success_flag, "target_gene_length": 0, "unrecovered_gene_length": 0,
                 "filter_reads_count": 0, "kmer_usage_rate": 0,
                 "seed": "", "assemble_kmer_size": _assemble_kmer_size_,
                 "kmer_limit": _kmer_limit_count_, "contig_coverage_depth": 0, "ref_avg_length": gene_avg_len}
    # read 产生的Kmer总数和不同种类Kmer的数量
    # ref_seq 产生的Kmer总数和不同种类Kmer的数量
    assemble_gene_info_dict[gene_name] = gene_info
    if _print_:
        print("Assembling the {} gene".format(gene_name))
    # 当read文件不存在时,返回False
    if not os.path.isfile(_reads_file_):
        print("Reads file of {} gene doesn't exist".format(gene_name))
        return assemble_gene_info_dict

    # 根据读长文件和参考序列文件 获取read_kmer_dict 和 ref_kmer_dict
    _read_kmer_dict_, _ref_kmer_dict_ = make_assemble_hash_dict(_reads_file_, _ref_file_, _assemble_kmer_size_,
                                                                _ref_reverse_complement_, True, _pos_, _print_=False)

    if not _read_kmer_dict_:
        # 当read_kmer_dict为0,认为是_short_contig_path_ 产生一个空文件
        with open(_short_contig_path_, 'w'):
            pass
        print("Reads file of {} gene doesn't have enough reads".format(gene_name))
        return assemble_gene_info_dict

    # # todo:修改为动态limit
    # if _kmer_limit_count_ < 0:
    #     _kmer_limit_count_ = util.dynamic_limit(_read_kmer_dict_, gene_avg_len, 2)

    # 根据kmer出现次数和是否在参考序列字典中 保留在参考序列字典中的kmer,删除出现数量低的kmer
    if _kmer_limit_count_ > 0:
        # kmer_dict中出现次数小于limit 且 在参考序列中没有出现 的kmer删除
        _filter_ = [_kmer_ for _kmer_ in _read_kmer_dict_ if _read_kmer_dict_[_kmer_][0] <= _kmer_limit_count_]
        for _kmer_ in _filter_:
            if _read_kmer_dict_[_kmer_][1] == 0:
                del _read_kmer_dict_[_kmer_]
        # 用于清理内存
        del _filter_

    # 获取reads数量信息
    filter_reads_count, _ = get_reads_info(_reads_file_)

    # 从ref_dict中获取种子 即同时出现在ref_dict和kmer_dict中的kmer
    # seed_list = [(x, _ref_kmer_dict_[x][0], _ref_kmer_dict_[x][1] / _ref_kmer_dict_[x][0]) for x in _ref_kmer_dict_
    #              if x in _read_kmer_dict_]

    # 设定种子的权重 选取在read中和ref中都存在的kmer
    # 种子的权重使用 ref_kmer的出现次数/ ref_kmer的出现的中位数 * read_kmer的出现次数/ read_kmer的出现的中位数
    shared_seed = [x for x in _ref_kmer_dict_ if x in _read_kmer_dict_]
    ref_kmer_count = [_ref_kmer_dict_[i][0] for i in shared_seed]
    read_kmer_count = [_read_kmer_dict_[i][0] for i in shared_seed]

    seed_list = [(x, _ref_kmer_dict_[x][0] / get_median(ref_kmer_count) * _read_kmer_dict_[x][0]
                  / get_median(read_kmer_count), _ref_kmer_dict_[x][1] / _ref_kmer_dict_[x][0]) for x in shared_seed]
    # 计算删除出现次数较少的kmer后剩余的Kmer总数
    read_kmer_count_sum = sum(_read_kmer_dict_[i][0] for i in _read_kmer_dict_)
    del shared_seed, ref_kmer_count, read_kmer_count

    # 根据kmer权重进行排序
    list.sort(seed_list, key=lambda x: x[1], reverse=True)

    # 当seed_dict为0,认为是_short_contig_path_ 产生一个空文件
    if not seed_list:
        with open(_short_contig_path_, 'w'):
            pass
        print("The ref file of {} gene does not have enough sequences".format(gene_name))
        return assemble_gene_info_dict

    # 选取出现次数最多的_kmer_作为起始种子
    # unique_kmer_number 计算不同的kmer的数量
    # used_kmer_count 使用的kmer的出现次数的总和
    _cur_seed_info_ = seed_list[0]
    contig, _, _, _, unique_kmer_number, used_kmer_count = \
        assemble_read_for_contig(_read_kmer_dict_, _cur_seed_info_[0], _iteration_)

    # 记录当前拼接出来最长的contig
    best_contig = contig
    best_kmer_count = used_kmer_count
    best_unique_kmer_number = unique_kmer_number
    best_seed_seq = _cur_seed_info_[0]

    short_contig_length = 0
    contig_length = 0
    # 如果序列太短 策略是更换seed
    if len(contig) / gene_avg_len < _min_percent_length_:
        # 记录当前contig和kmer相关的信息
        change_count = 0
        last_seed_info = _cur_seed_info_
        # 从seed_list中选取一个seed的信息元组
        for new_seed_info in seed_list[1:]:
            # 不选择和当前seed有1-2bp重合的kmer
            if last_seed_info[0][1:] == new_seed_info[0][:-1] or last_seed_info[0][2:] == new_seed_info[0][:-2]:
                last_seed_info = new_seed_info
                continue
            # 记录seed信息
            last_seed_info = new_seed_info
            # print("The contig length is {} and is shorter than avg_len. Use another seed to re-assemble."
            #       .format(len(contig)))
            # 产生的新的contig
            contig, _, _, _, unique_kmer_number, used_kmer_count = \
                assemble_read_for_contig(_read_kmer_dict_, new_seed_info[0], _iteration_)
            change_count += 1
            # 如果新的contig的长度比之前存储的best_contig长,则更改记录的信息
            if len(contig) > len(best_contig):
                best_contig, best_seed_seq, best_unique_kmer_number, best_kmer_count = \
                    contig, new_seed_info[0], unique_kmer_number, used_kmer_count
            if len(best_contig) / gene_avg_len >= _min_percent_length_ or change_count >= int(_change_seed_):
                break

        # 当最后获得的best_contig的长度较短:将文件存储进入_short_contig_path_ 否则存入_contig_path_
        # 当获得的best_contig长度达标,则is_normal_contig:True assemble_success_flag:True
        if len(best_contig) / gene_avg_len >= _min_percent_length_:
            contig_length = len(best_contig)
            is_normal_contig = True
            assemble_success_flag = True
        else:
            short_contig_length = len(best_contig)
    # 如果序列太长,策略是增大_kmer_limit_count_
    elif len(contig) / gene_avg_len > _max_percent_length_:
        is_normal_contig = True
        assemble_success_flag = True
        # 加大limit
        for temp_kmer_limit in range(_kmer_limit_count_ + 1, 32):
            # 对filtered_dict再次精简
            _filter_ = [x for x in _read_kmer_dict_ if _read_kmer_dict_[x][0] <= temp_kmer_limit]
            for x in _filter_:
                if _read_kmer_dict_[x][1] == 0:
                    del _read_kmer_dict_[x]
            del _filter_
            contig, _, _, _, unique_kmer_number, used_kmer_count \
                = assemble_read_for_contig(_read_kmer_dict_, _cur_seed_info_[0])
            if _max_percent_length_ >= len(contig) / gene_avg_len >= _min_percent_length_:
                best_contig, best_unique_kmer_number, best_kmer_count = contig, unique_kmer_number, used_kmer_count
                break
            # 当增大_kmer_limit_count_后,若是序列低于标准,则选用当前的best_contig
            elif len(contig) / gene_avg_len < _min_percent_length_:
                break
            # 当增大_kmer_limit_count_后,序列仍较长时,则更换当前最优contig
            else:
                best_contig, best_unique_kmer_number, best_kmer_count = contig, unique_kmer_number, used_kmer_count
        contig_length = len(best_contig)
    # 序列长度合理,则直接选用当前的best_contig
    else:
        contig_length = len(best_contig)
        is_normal_contig = True
        assemble_success_flag = True

    if is_normal_contig:
        if _print_:
            print("Assemble {} gene succeed. Ref length: {}, best contig length: {}.".format
                  (gene_name, gene_avg_len, len(best_contig)))
        with open(_contig_path_, 'w') as out:
            out.write(
                '>' + os.path.split(_out_dir_)[-1] + "_" + gene_name + '_recovered_k' + str(_assemble_kmer_size_) +
                "_" + str(len(best_contig)) + '\n')
            out.write(best_contig + '\n')
    else:
        if _print_:
            print("Assemble {} gene failed. Ref length: {}, best contig length: {}.".format
                  (gene_name, gene_avg_len, len(best_contig)))
        with open(_short_contig_path_, 'w') as out:
            out.write('>' + os.path.split(_out_dir_)[-1] + "_" + gene_name + '_unrecovered_k' +
                      str(_assemble_kmer_size_) + "_" + str(len(best_contig)) + '\n')
            out.write(best_contig + '\n')

    # 使用used_kmer_count*kmer_size / unique_kmer_count 计算contig的覆盖率
    # 覆盖深度
    contig_coverage_depth = round(best_kmer_count * _assemble_kmer_size_ / best_unique_kmer_number, 2)

    gene_info["assemble_success_flag"] = assemble_success_flag
    gene_info["target_gene_length"] = contig_length
    gene_info["unrecovered_gene_length"] = short_contig_length
    gene_info["filter_reads_count"] = filter_reads_count
    gene_info["kmer_usage_rate"] = round(best_kmer_count / read_kmer_count_sum, 2)
    gene_info["seed"] = best_seed_seq
    gene_info["contig_coverage_depth"] = contig_coverage_depth
    assemble_gene_info_dict[gene_name] = gene_info
    return assemble_gene_info_dict


def assemble_flow(_input_read_path_: str, _out_dir_: str, _ref_path_: str, _assemble_kmer_size_: int,
                  _assemble_thread_: int = 4, _ref_reverse_complement_: bool = True, _pos_: bool = True,
                  _change_seed_: int = 1000, _kmer_limit_count_: int = 2, _min_percent_length_: float = 1.0,
                  _max_percent_length_: float = 2.0, _iteration_: int = 1000):
    # 设定assemble contig的输出文件夹
    # 当存在这些目录文件时,会删除该目录及其下的所有文件,而后重新生成文件夹
    assemble_out_dir = os.path.join(_out_dir_, "target_genes")
    if not os.path.isdir(assemble_out_dir):
        os.makedirs(assemble_out_dir)
    else:
        rmtree(assemble_out_dir)
        os.makedirs(assemble_out_dir)

    if not os.path.isdir(os.path.join(assemble_out_dir, "unrecovered_genes")):
        os.makedirs(os.path.join(assemble_out_dir, "unrecovered_genes"))

    print('Assemble reads')
    file_dict = collections.defaultdict(dict)
    for ref_file_path in get_file_list(_ref_path_):
        gene_name = os.path.basename(ref_file_path).split(".")[0]
        file_dict[gene_name]["ref_file"] = ref_file_path
        if os.path.isdir(_input_read_path_) and os.path.isfile(os.path.join(_input_read_path_, gene_name + ".fasta")):
            file_dict[gene_name]["reads_file"] = os.path.join(_input_read_path_, gene_name + ".fasta")
        elif os.path.isfile(_input_read_path_) and gene_name in _input_read_path_:
            file_dict[gene_name]["reads_file"] = _input_read_path_
        else:
            file_dict[gene_name]["reads_file"] = None
    # 用于存储基因是否过滤成功的信息
    assemble_gene_info_dict = collections.defaultdict(dict)
    count = 0
    with ThreadPoolExecutor(max_workers=min(_assemble_thread_, len(file_dict))) as executor:
        futures = [executor.submit(
            partial(assemble_one_file, _out_dir_=assemble_out_dir, _assemble_kmer_size_=_assemble_kmer_size_,
                    _ref_reverse_complement_=_ref_reverse_complement_, _pos_=_pos_,
                    _change_seed_=_change_seed_, _kmer_limit_count_=_kmer_limit_count_,
                    _min_percent_length_=_min_percent_length_, _max_percent_length_=_max_percent_length_,
                    _iteration_=_iteration_),
            _reads_file_=file_info["reads_file"], _ref_file_=file_info["ref_file"])
            for gene_name, file_info in file_dict.items()]
        for future in as_completed(futures):
            assemble_gene_info_dict.update(future.result())
            count += 1
            if count % 10 == 0:
                print("INFO: {} / {} genes have been assembled!".format(count, len(file_dict)))

    assemble_flag_list = [i["assemble_success_flag"] for i in assemble_gene_info_dict.values()]
    print("Assemble Done! {} / {} succeed".format(assemble_flag_list.count(True), len(assemble_flag_list)))
    # 输出log信息
    log_file = os.path.join(assemble_out_dir, "assemble_log.csv")
    log(log_file, "gene_id", "ref_avg_length", "target_gene_length", "unrecovered_gene_length", "contig_coverage_depth",
        "assemble_kmer", "assemble_seed", "kmer_limit", "filter_reads_count",
        "kmer_usage_rate")
    for gene_name, gene_info in assemble_gene_info_dict.items():
        log(log_file, gene_name, gene_info["ref_avg_length"], gene_info["target_gene_length"],
            gene_info["unrecovered_gene_length"],
            gene_info["contig_coverage_depth"], gene_info["assemble_kmer_size"],
            gene_info["seed"], gene_info["kmer_limit"], gene_info["filter_reads_count"],
            gene_info["kmer_usage_rate"])
    return assemble_gene_info_dict


if __name__ == "__main__":
    # todo: 评估与参考序列的相似性 用参与拼接的K-mer在参考序列中出现的平均次数替代
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Miner
    zzhen@sculab ''')
    pars.add_argument('-i', dest='input_dir', help='Input a file(directory) with reads data', required=True)
    pars.add_argument("-r", dest="reference", type=str, help="Input a file(directory) with references", required=True)
    pars.add_argument("-o", dest="output_dir", type=str, help="Output directory.", required=False,
                      default="easy353_output")
    pars.add_argument("-k", dest="assemble_kmer", type=int, help="Kmer setting for assembling reads. Default:41",
                      default=41)
    pars.add_argument("-s", dest="step_length", type=int,
                      help="Step length of the sliding window on the reads. Default:1", default=1)
    pars.add_argument("-t", dest="assemble_thread", type=int,
                      help="Threads setting for assembling reads. Default:4", default=4)
    pars.add_argument("-kmer_limit", dest="kmer_limit", type=int, help="Limit of kmer count. Default:3", default=2)
    pars.add_argument("-min", dest="minimum_length_ratio", type=float,
                      help="The minimum ratio of contig length to reference average length. Default:1.0", default=1.0)
    pars.add_argument("-max", dest="maximum_length_ratio", type=float,
                      help="The maximum ratio of contig length to reference average length. Default:2.0", default=2.0)
    pars.add_argument("-change_seed", dest="change_seed", type=int, help="Times of changing seed. Default:32",
                      default=32)
    pars.add_argument("-fast", dest="fast", action="store_true", help="Whether to use fast mode.")
    # 默认限制使用的参考序列的数量
    pars.add_argument("-reference_number", dest="reference_number", type=int,
                      help="The number of the reference sequences used to build hash table. Default:all", default=None)
    args = pars.parse_args()

    assemble_flow(_input_read_path_=args.input_dir, _out_dir_=args.output_dir, _ref_path_=args.reference,
                  _assemble_kmer_size_=args.assemble_kmer,
                  _assemble_thread_=args.assemble_thread, _ref_reverse_complement_=True, _pos_=True,
                  _change_seed_=args.change_seed, _kmer_limit_count_=args.kmer_limit,
                  _min_percent_length_=args.minimum_length_ratio, _max_percent_length_=args.maximum_length_ratio,
                  _iteration_=1000)
    # assemble_one_file("/Users/zzhen/Desktop/test/7367.fasta", "/Users/zzhen/Desktop/debug"
    #                   , "/Users/zzhen/Desktop/new/353gene/7367.fasta",
    #                   41, True, True, 32, 2, 1.0, 2.0, 1000, False)
