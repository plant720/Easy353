#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:50 下午
# @Author     : zzhen
# @File       : interface.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import gc
import os
import re
import shutil
import time
import psutil
import threading
import collections
from urllib.error import HTTPError
from urllib.request import urlretrieve
import PySimpleGUI as sg
from multiprocessing import Pool
from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed
# import my python file
from src.filter import do_reads_filter, do_pair_reads_filter, refilter_one_gene
from src.utils import make_ref_kmer_dict, get_file_size, get_ref_info, \
    get_reads_info, log, get_file_list
from src.build_database import generate_download_info, detect_classification
from src.build_database import generate_fasta_path, generate_gene_file
from src.assemble import assemble_one_file


# 进行窗口的的初始化
def window_init():
    sg.theme('Default 1')
    # Parameters
    # basic_para_frame
    basic_para_frame = [
        # sequence_data
        [sg.Text('Single fq file', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-q-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with unpaired (single-end) reads.",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-q-')],

        # paired-end data
        [sg.Text('Forward fq file', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with forward paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-1-')],
        # paired-end data
        [sg.Text('Reverse fq file', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with reverse paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-2-')],
        # classification
        [sg.Text('Classification', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key="-classification-", size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input the interested classification.", readonly=False, background_color='#BAD9BC')],
        # output
        [sg.Text('Output dir', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-o-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip='Output directory.',
                  readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Times", 12), target='-o-')],
    ]
    # general_para_frame
    general_para_frame = [
        [sg.Text('Exclude file', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-exclude_file-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip='The file documenting species that need to be excluded',
                  readonly=False, background_color='#BFE9FF'),
         sg.FileBrowse("File", font=("Times", 12), target='-exclude_file-')],

        [sg.Text("Exclude", size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(key="-exclude-", size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input species that need to be excluded", readonly=False,
                  background_color='#BFE9FF')],

        # kmer for filter and assembly
        [sg.Text('K-mer for filter', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(31, key='-k1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Kmer setting for filtering reads. [Default:31]"),
         sg.Text('K-mer for assembly', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(41, key='-k2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Kmer setting for assembling reads. [Default:41]"
                  )],

        # threads for filter and assembly
        [sg.Text('Threads for filtering', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(1, key='-t1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Threads setting for filtering reads. [Default:4]"),
         sg.Text('Threads for assembly', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(4, key='-t2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Threads setting for assembling reads. [Default:4]"
                  )],
    ]
    advanced_para_frame = [
        # sliding_window
        [sg.Text('Step length', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(1, key='-s-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Step length of the sliding window on the reads. [Default:1]"
                  ),
         sg.Text('Ref number', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(100, key='-reference_number-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The number of the reference sequences used to build hash table. [Default:100]")],
        # change_seed and limit_count
        [sg.Text('Change seed', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(32, key='-change_seed-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Times of changing seed. [Default:32]"),

         sg.Text('Kmer limit', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(2, key='-kmer_limit-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Limit of kmer count. [Default:2]"
                  )],
        # min_percent_length and max_percent_length
        [sg.Text('Min length ratio', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(1.0, key='-minimum_length_ratio-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The minimum ratio of contig length to reference average length. [Default:1.0]"),

         sg.Text('Max length ratio', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(2.0, key='-maximum_length_ratio-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The maximum ratio of contig length to reference average length. [Default:2.0]")
         ],

        # reverse_complement
        [sg.Text('Ref reverse complement', size=(21, 1), justification='center', font=("Times", 12)),
         sg.Checkbox(text="True", default=True,
                     key='-reverse_complement-', size=(12, 1), font=("Times", 12), ),

         sg.Text('Generate scaffold', size=(15, 1), justification='center', font=("Times", 12)),
         sg.Checkbox(text="True", default=False, key='-generate_scaffold-', size=(12, 1), font=("Times", 12), )
         ]
    ]

    # left_col
    left_col = [
        # 基本设定
        [sg.Text('Parameters', justification='center', font=("Times", 14, 'bold'))],
        [sg.Frame('Basic', layout=basic_para_frame, expand_x=True, font=("Times", 13, 'bold'))],
        [sg.Frame('General', layout=general_para_frame, expand_x=True, font=("Times", 13, 'bold'))],
        [sg.Frame('Advanced', layout=advanced_para_frame, expand_x=True, font=("Times", 13, 'bold'))],
    ]
    # right_col
    right_col = [
        # [sg.Output(size=(100, 25), font=("Times", 12, 'bold'), text_color="black", key='-OUTPUT-')],
        [sg.Text('Running log', justification='center', font=("Times", 14, 'bold'))],
        [sg.Multiline(size=(60, 20), font=("Times", 12), expand_x=True, expand_y=True,
                      key='-OUTPUT-', autoscroll=True, reroute_stdout=True, reroute_stderr=True,
                      auto_refresh=True)],
        # 按钮
        [sg.Button('Run', font=("Times", 12), auto_size_button=True),
         sg.Button('Reset', font=("Times", 12), auto_size_button=True),
         sg.Button('Close', font=("Times", 12), auto_size_button=True)],
    ]

    info_frame = [
        [sg.Text('Easy353', justification="center", font=("Times", 24, "bold"))],
        [sg.Text('Authors: zzhen', font=("Times", 12)),
         sg.Text('Email: zzhen0302@163.com', font=("Times", 12))],
        [sg.Text('Version: 1.2.0', font=("Times", 12))]
    ]
    # 设定最终布局
    layout = [
        [sg.Push(), sg.Column(info_frame, element_justification='c'), sg.Push()],
        [sg.Pane([sg.Column(left_col, element_justification='c', expand_x=True, expand_y=True),
                  sg.Column(right_col, element_justification='c', expand_x=True, expand_y=True)],
                 orientation='horizontal', relief=sg.RELIEF_SUNKEN, key='-PANE-')],
    ]

    window = sg.Window(title='Easy353', layout=layout)
    return window


# 用于下载文件的函数
def download_files(_spec_info_: dict, output_dir: str):
    info = None
    url = _spec_info_['Fasta file url']
    file_path = os.path.join(output_dir, _spec_info_['Fasta file name'])
    if os.path.isfile(file_path):
        info = "INFO: File {} already exists".format(file_path)
    else:
        try:
            urlretrieve(url, file_path)
        except HTTPError:
            info = "INFO: Url {} does not exist".format(url)
    return info


# 用于根据传入的分类名称和输出文件夹确定需要下载的文件
def download_species_thread(output_dir: str, classification: list, window: sg.Window, max_threads: int = 8):
    print("INFO: Start of generating reference...")
    if output_dir is not None:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    # 当需要下载文件时
    print("INFO: Get species information that needs to be downloaded")
    # 根据传入的分类信息需要下载的物种信息 是元素为dict的list
    down_spec_info = generate_download_info(detect_classification(classification))
    # 检测已经存在的文件名
    existed_files = [os.path.basename(i) for i in generate_fasta_path(output_dir)]
    # 获取还未下载的文件
    down_spec_info = [i for i in down_spec_info if i["Fasta file name"] not in existed_files]
    print("INFO: Download species data")
    print("INFO: Total {} species need to be downloaded".format(len(down_spec_info)))
    # 下载文件
    count = 1
    # 使用concurrent.features中的ThreadPoolExecutor类来实现多线程下载
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
        # 可以使用偏函数构造出单参数函数 一般来说必须采用关键字参数的形式给构建出的偏函数传参，因为如果以无关键字参数的方式，该实参将试图传递给第一个形参
        # http://c.biancheng.net/view/5674.html
        futures = [executor.submit(partial(download_files, output_dir=output_dir), _spec_info_) for _spec_info_ in
                   down_spec_info]
        for future in as_completed(futures):
            if future.result() is not None:
                print(future.result())
            if count % 10 == 0:
                print("INFO: {} / {} has been downloaded".format(count, len(down_spec_info)))
            count += 1
    print("INFO: Download finished!")
    # 向主进程进行通信
    window.write_event_value('-DOWNLOAD DONE-', "-DOWNLOAD DONE-")


# 根据下载的物种353数据和需要排除的物种信息 生成参考序列
def generate_ref_thread(output_dir: str, exclude: list, exclude_file: str, window: sg.Window):
    # 根据已经下载的文件夹
    exclude_species = []
    if exclude is not None:
        exclude_species.extend(exclude)
    if exclude_file is not None:
        if not os.path.isfile(exclude_file):
            print("INFO: The exclude file does not exist!")
        else:
            with open(exclude_file, "r") as _file_:
                _exclude_ = _file_.readlines()
                _exclude_ = [x.strip() for x in _exclude_ if x.strip()]
                exclude_species.extend(_exclude_)
    exclude_species = [species.capitalize().replace(" ", "_") for species in exclude_species]
    print("INFO: Generating species data into reference files")
    file_path_list = generate_fasta_path(output_dir)
    generate_gene_file(_file_path_list_=file_path_list,
                       _output_dir_=os.path.join(output_dir, "353gene"), _exclude_species_=exclude_species)
    print("INFO: End of generating reference!")
    window.write_event_value("-GENERATE DONE-", "-GENERATE DONE-")


# 根据参考序列文件 将测序数据中与参考序列相关的read过滤出来
def filter_read_thread(window: sg.Window, _read_data_tuple_: tuple, _out_dir_: str, _reference_path_: str,
                       _kmer_size_: int, _step_size_: int = 1, _ref_reverse_complement_: bool = False,
                       _read_reverse_complement_: bool = False, refilter: bool = True,
                       _clear_: bool = True, _pos_: bool = True, _paired_reads_: bool = True,
                       _thread_for_filter_: int = 4, _print_: bool = True, _ref_number_: int = None):
    print("INFO: Start of reads filtering")
    # 初始化操作
    if not os.path.isdir(_out_dir_):
        os.makedirs(_out_dir_)
    # 当参考序列不进行反向互补时，read进反向互补
    if not _ref_reverse_complement_:
        _read_reverse_complement_ = True
    _time_build_hash_start_ = time.time()
    # ref_length_dict 记录reference文件的read平均长度
    # ref_path_dict 记录reference文件的路径信息
    ref_length_dict, ref_path_dict = get_ref_info(_reference_path_)
    # 构建哈西字典
    print("INFO: Building hash table")
    # 处理reference 并生成具体kmer存储在ref_kmer_dict
    # 用于记录记录所有的reference信息  字典的key是kmer seq value是一个list,
    # list[0] 是kmer出现次数
    # list[1] 是kmer出现在某条参考基因上的百分比位置之和
    # list[2:]是该kmer在哪些文件(基因)中出现
    print("INFO: The ref number used for one gene is {}".format(_ref_number_))
    ref_kmer_dict = make_ref_kmer_dict(_reference_path_, _kmer_size_, _ref_reverse_complement_, _pos_, _print_,
                                       _ref_number_)
    print("INFO: Hash table of K-mers has been made. Time used: {} s".format(
        round(time.time() - _time_build_hash_start_), 2))

    # _clear_用于将输出文件夹中的原有同名文件删除,新建文件
    if _clear_:
        for key in ref_path_dict:
            with open(os.path.join(_out_dir_, key + ".fasta"), 'w'):
                pass
    print("INFO: The memory usage of the current process is: {} GB".
          format(round(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, 2)))
    print('INFO: Filter reads from fq_files based on hash table with {} thread'.format(_thread_for_filter_))
    _time_filter_start_ = time.time()
    # 开始进行read过滤
    # 当测序数据是unpaired_end时
    if not _paired_reads_:
        _unpaired_file_ = _read_data_tuple_[0]
        # 当建立的hash表的内存超过系统内存一半时，使用单线程进行过滤
        do_reads_filter(ref_kmer_dict, _kmer_size_, _step_size_, _unpaired_file_,
                        _out_dir_, _read_reverse_complement_, _print_)
        # count = 1
        # # 使用concurrent.features中的ThreadPoolExecutor类来实现多线程下载
        # with ThreadPoolExecutor(max_workers=_thread_for_filter_) as executor:
        #     # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
        #     # 可以使用偏函数构造出单参数函数
        #     futures = [executor.submit(partial(do_reads_filter, _ref_kmer_dict_=ref_kmer_dict,
        #                                        _kmer_size_=_kmer_size_, _step_size_=_step_size_,
        #                                        _out_dir_=_out_dir_,
        #                                        _read_reverse_complement_=_read_reverse_complement_,
        #                                        _print_=_print_), unpaired_file)
        #                for unpaired_file in _unpaired_file_list_]
        #     for future in as_completed(futures):
        #         print("INFO: {} / {} unpaired files has been downloaded".format(count, len(_unpaired_file_list_)))
        #         count += 1
    # 当测序数据是paired-end时
    if _paired_reads_:
        paired_file_forward = _read_data_tuple_[0]
        paired_file_reverse = _read_data_tuple_[1]
        # 处理根据内存占比，选用多进程或是单进程
        if _thread_for_filter_ == 1:
            # 当设定为单线程
            do_pair_reads_filter(ref_kmer_dict, _kmer_size_, _step_size_,
                                 paired_file_forward, paired_file_reverse,
                                 _out_dir_, 0, 1, _read_reverse_complement_, _print_)
        else:
            # 创建进程池
            print("INFO: This step will take about 20 minutes and no log information is displayed")
            pool = Pool(_thread_for_filter_)
            # 创建进程池中的进程
            # 多进程读取一个文件
            for i in range(_thread_for_filter_):
                pool.apply_async(func=do_pair_reads_filter,
                                 args=(ref_kmer_dict, _kmer_size_, _step_size_,
                                       paired_file_forward, paired_file_reverse,
                                       _out_dir_, i, _thread_for_filter_, _read_reverse_complement_, _print_,))
            pool.close()
            pool.join()
            # # 使用concurrent.features中的ThreadPoolExecutor类来实现多线程下载
            # with ThreadPoolExecutor(max_workers=_thread_for_filter_) as executor:
            #     # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
            #     # 可以使用偏函数构造出单参数函数
            #     futures = [executor.submit(partial(do_pair_reads_filter, _ref_kmer_dict_=ref_kmer_dict,
            #                                        _kmer_size_=_kmer_size_, _step_size_=_step_size_,
            #                                        fq_file_1=paired_file_forward, fq_file_2=paired_file_reverse,
            #                                        _out_dir_=_out_dir_, t_count=_thread_for_filter_,
            #                                        _read_reverse_complement_=_read_reverse_complement_,
            #                                        _print_=_print_), t_id=i) for i in range(_thread_for_filter_)]
    print("INFO: Time used for filtering reads: {} s".format(round(time.time() - _time_filter_start_), 2))
    # 过滤完后 手动清理内存
    del ref_kmer_dict
    gc.collect()

    # 根据过滤后的数据 获取过滤出的reads数量和reads总长度
    filter_gene_info_dict = collections.defaultdict(dict)
    for gene_name, ref_avg_length in ref_length_dict.items():
        filter_gene_info_dict[gene_name]["ref_avg_length"] = ref_avg_length
        filter_gene_info_dict[gene_name]["filter_kmer"] = _kmer_size_
        filter_gene_info_dict[gene_name]["filter_reads_count"], filter_gene_info_dict[gene_name][
            "filter_reads_length"] = get_reads_info(os.path.join(_out_dir_, gene_name + ".fasta"))

    # 当需要re-filter时
    if refilter:
        # 记录需要过滤的基因 对应参考基因平均长度和参考序列所在文件 key:gene_name value:[avg_len,gene_path]
        refilter_gene_info_dict = collections.defaultdict(dict)
        if not os.path.isdir(os.path.join(_out_dir_, 'big_reads')):
            os.makedirs(os.path.join(_out_dir_, 'big_reads'))
        # ref_length_dict 记录reference文件的read平均长度
        for gene_name, gene_info in filter_gene_info_dict.items():
            if os.path.isfile(os.path.join(_out_dir_, gene_name + ".fasta")):
                # 如果覆盖深度>512x 或过滤后的文件>8m
                # 该处覆盖深度采用简单的 c = LN / G
                if (gene_info["filter_reads_length"] / gene_info.get("ref_avg_length") > 512) \
                        and get_file_size(os.path.join(_out_dir_, gene_name + ".fasta")) > 8:
                    # 将需要过滤的文件移动到一个文件夹中
                    shutil.move(os.path.join(_out_dir_, gene_name + ".fasta"),
                                os.path.join(_out_dir_, 'big_reads', gene_name + ".fasta"))
                    refilter_gene_info_dict[gene_name] = gene_info
        # 再次清理内存
        gc.collect()
        if not refilter_gene_info_dict:
            print("INFO: Big reads files don't exist, no need to re-filter!")
        else:
            print("INFO: re-filter")
            # 使用多线程进行过滤
            with ThreadPoolExecutor(max_workers=min(_thread_for_filter_, len(refilter_gene_info_dict))) as executor:
                futures = [executor.submit(
                    partial(refilter_one_gene, _refilter_kmer_size_=_kmer_size_ + 2, _out_dir_=_out_dir_,
                            _refilter_step_size_=1), gene_name=gene_name, gene_avg_len=gene_info["ref_avg_length"],
                    _reference_path_=ref_path_dict[gene_name])
                    for gene_name, gene_info in refilter_gene_info_dict.items()]
                for future in as_completed(futures):
                    refilter_gene_info_dict.update(future.result())
                    # 输出log信息
        filter_gene_info_dict.update(refilter_gene_info_dict)
        log_file = os.path.join(_out_dir_, "filter_log.csv")
        log(log_file, "gene_id", "ref_avg_length", "filter_reads_count", "filter_kmer")
        for gene_name, gene_info in filter_gene_info_dict.items():
            log(log_file, gene_name, gene_info["ref_avg_length"], gene_info["filter_reads_count"],
                gene_info["filter_kmer"])
        print("INFO: End of reads filtering!")
        window.write_event_value("-FILTER DONE-", "-FILTER DONE-")


def assemble_read_thread(window: sg.Window, _input_read_path_: str, _out_dir_: str, _ref_path_: str,
                         _assemble_kmer_size_: int, _assemble_thread_: int = 4, _ref_reverse_complement_: bool = True,
                         _pos_: bool = True, _change_seed_: int = 1000, _kmer_limit_count_: int = 2,
                         _min_percent_length_: float = 1.0, _max_percent_length_: float = 2.0,
                         _iteration_: int = 1000, _write_scaffold_: bool = True):
    print('INFO: Start of reads assembly')
    # 设定assemble contig的输出文件夹
    # 当存在这些目录文件时,会删除该目录及其下的所有文件,而后重新生成文件夹
    if not os.path.isdir(os.path.join(_out_dir_, 'contigs')):
        os.makedirs(os.path.join(_out_dir_, 'contigs'))
    else:
        shutil.rmtree(os.path.join(_out_dir_, 'contigs'))
        os.makedirs(os.path.join(_out_dir_, 'contigs'))
    if not os.path.isdir(os.path.join(_out_dir_, 'short_contigs')):
        os.makedirs(os.path.join(_out_dir_, 'short_contigs'))
    else:
        shutil.rmtree(os.path.join(_out_dir_, 'short_contigs'))
        os.makedirs(os.path.join(_out_dir_, 'short_contigs'))
    if not os.path.isdir(os.path.join(_out_dir_, 'scaffolds')):
        os.makedirs(os.path.join(_out_dir_, 'scaffolds'))
    else:
        shutil.rmtree(os.path.join(_out_dir_, 'scaffolds'))
        os.makedirs(os.path.join(_out_dir_, 'scaffolds'))

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
    # 清除没有reads的文件
    for gene_name, file_info in file_dict.items():
        if file_info["reads_file"] is None:
            file_dict.pop(gene_name)

    count = 0
    with ThreadPoolExecutor(max_workers=min(_assemble_thread_, len(file_dict))) as executor:
        futures = [executor.submit(
            partial(assemble_one_file, _out_dir_=_out_dir_, _assemble_kmer_size_=_assemble_kmer_size_,
                    _ref_reverse_complement_=_ref_reverse_complement_, _pos_=_pos_,
                    _change_seed_=_change_seed_, _kmer_limit_count_=_kmer_limit_count_,
                    _min_percent_length_=_min_percent_length_, _max_percent_length_=_max_percent_length_,
                    _iteration_=_iteration_, _generate_scaffold_=_write_scaffold_, _print_=False),
            _reads_file_=file_info["reads_file"], _ref_file_=file_info["ref_file"])
            for gene_name, file_info in file_dict.items()]
        for future in as_completed(futures):
            assemble_gene_info_dict.update(future.result())
            count += 1
            if count % 10 == 0:
                print("INFO: {} / {} genes have been assembled!".format(count, len(file_dict)))

    assemble_flag_list = [i["assemble_success_flag"] for i in assemble_gene_info_dict.values()]
    print("INFO: Assemble Done! {} / {} succeed".format(assemble_flag_list.count(True), len(assemble_flag_list)))
    # 输出log信息
    log_file = os.path.join(_out_dir_, "assemble_log.csv")
    log(log_file, "gene_id", "contig_length", "short_contig_length", "scaffold_length", "contig_coverage_depth",
        "assemble_kmer", "assemble_seed", "kmer_limit", "filter_reads_count",
        "kmer_usage_rate")
    for gene_name, gene_info in assemble_gene_info_dict.items():
        log(log_file, gene_name, gene_info["contig_length"], gene_info["short_contig_length"],
            gene_info["scaffold_length"], gene_info["contig_coverage_depth"], gene_info["assemble_kmer_size"],
            gene_info["seed"], gene_info["kmer_limit"], gene_info["filter_reads_count"],
            gene_info["kmer_usage_rate"])
    print("INFO: End of reads assembly!")
    window.write_event_value("-ASSEMBLE DONE-", "-ASSEMBLE DONE-")


def easy353_gui():
    window = window_init()
    fastq_files, _paired_reads_ = tuple(), False
    fq_file_1, fq_file_2, unpaired_fq_file, output_dir = None, None, None, None
    filter_kmer, assemble_kmer, filter_thread, assemble_thread = 31, 41, 1, 4
    step_length, kmer_limit, minimum_length_ratio, maximum_length_ratio, change_seed, reference_number \
        = 1, 2, 1.0, 2.0, 32, 100
    # fast 即 ref_reverse_complement
    generate_scaffold, fast = False, False
    # get_ref parameter
    classification, exclude_file, exclude = None, None, None

    thread = None
    while True:
        # filter and assemble parameter
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Close':  # if user closes window or clicks cancel
            break
        if event == 'Reset' and not thread:
            window['-1-'].update('')
            window['-2-'].update('')
            window['-q-'].update('')
            window['-o-'].update('')
            window['-classification-'].update('')

            window['-exclude_file-'].update('')
            window['-exclude-'].update('')
            window['-k1-'].update(31)
            window['-k2-'].update(41)
            window['-t1-'].update(1)
            window['-t2-'].update(4)
            window['-s-'].update(1)
            window['-reference_number-'].update(100)

            window['-change_seed-'].update(32)
            window['-kmer_limit-'].update(2)
            window['-minimum_length_ratio-'].update(1.0)
            window['-maximum_length_ratio-'].update(2.0)
            window['-generate_scaffold-'].update(False)
            window['-reverse_complement-'].update(True)

            fastq_files, _paired_reads_ = tuple(), False
            fq_file_1, fq_file_2, unpaired_fq_file, output_dir = None, None, None, None
            filter_kmer, assemble_kmer, filter_thread, assemble_thread = 31, 41, 1, 1
            step_length, kmer_limit, minimum_length_ratio, maximum_length_ratio, change_seed, reference_number \
                = 1, 2, 1.0, 2.0, 32, 100
            # fast 即 ref_reverse_complement
            generate_scaffold, fast = False, False
            # get_ref parameter
            classification, exclude_file, exclude = None, None, None

        if event == 'Run' and not thread:
            # 基本参数部分
            if values["-1-"]:
                fq_file_1 = values["-1-"]
            if values["-2-"]:
                fq_file_2 = values["-2-"]
            if values["-q-"]:
                unpaired_fq_file = values["-q-"]
            if unpaired_fq_file or (fq_file_1 and fq_file_2):
                if unpaired_fq_file:
                    fastq_files = (unpaired_fq_file,)
                if fq_file_1 and fq_file_2:
                    fastq_files = (fq_file_1, fq_file_2)
                    _paired_reads_ = True
            else:
                sg.Popup("Please input fastq file(s)!",
                         title='Info', keep_on_top=True, font=("Times", 12))
                continue
            if values["-o-"]:
                output_dir = values["-o-"]
            else:
                sg.Popup("Please input an output directory!",
                         title='Info', keep_on_top=True, font=("Times", 12))
                continue
            if values["-classification-"]:
                # 处理classification
                classification = values["-classification-"]
                classification = re.split(r"[,，、 ]", classification)
            else:
                sg.Popup("Please input a classification!", title='Info', keep_on_top=True, font=("Times", 12))
                continue

            if values["-exclude-"]:
                exclude = values["-exclude-"]
                exclude = re.split(r"[,，、 ]", exclude)
            if values["-exclude_file-"]:
                exclude_file = values["-exclude_file-"]

            if values["-k1-"]:
                filter_kmer = int(values["-k1-"])
            if values["-k2-"]:
                assemble_kmer = int(values["-k2-"])
            if values["-t1-"]:
                filter_thread = int(values["-t1-"])
            if values["-t2-"]:
                assemble_thread = int(values["-t2-"])
            if values["-s-"]:
                step_length = int(values["-s-"])
            if values["-reference_number-"]:
                reference_number = int(values["-reference_number-"])
            if values["-change_seed-"]:
                change_seed = int(values["-change_seed-"])
            if values["-kmer_limit-"]:
                kmer_limit = int(values["-kmer_limit-"])
            if values["-minimum_length_ratio-"]:
                minimum_length_ratio = float(values["-minimum_length_ratio-"])
            if values["-maximum_length_ratio-"]:
                maximum_length_ratio = float(values["-maximum_length_ratio-"])
            if values["-reverse_complement-"]:
                fast = True
            if values["-generate_scaffold-"]:
                generate_scaffold = True

            # 需要处理参考序列文件
            # 进入物种数据文件下载的进程
            thread = threading.Thread(
                target=download_species_thread,
                args=(output_dir, classification, window),
                daemon=True)
            thread.start()

        if event == '-DOWNLOAD DONE-':
            # 当数据下载成功后，进入建立参考序列的新进程
            # thread变量发生变化
            thread = threading.Thread(
                target=generate_ref_thread,
                args=(output_dir, exclude, exclude_file, window),
                daemon=True)
            thread.start()

        if event == '-GENERATE DONE-':
            # 当参考序列获取成功后，进入读长过滤的新进程
            reference = os.path.join(output_dir, "353gene")
            filter_out = os.path.join(output_dir, "easy353", "reads")
            if not os.path.isdir(filter_out):
                os.makedirs(filter_out)
            # 当数据下载成功后，进入建立参考序列的新进程
            # thread变量发生变化
            thread = threading.Thread(
                target=filter_read_thread,
                args=(window, fastq_files, filter_out, reference, filter_kmer, step_length, fast,
                      False, True, True, True, _paired_reads_, filter_thread, False,
                      reference_number),
                daemon=True)
            thread.start()

        if event == "-FILTER DONE-":
            reference = os.path.join(output_dir, "353gene")
            filter_out = os.path.join(output_dir, "easy353", "reads")
            assemble_out = os.path.join(output_dir, "easy353")
            if not os.path.isdir(assemble_out):
                os.makedirs(assemble_out)
            thread = threading.Thread(
                target=assemble_read_thread,
                args=(window, filter_out, assemble_out, reference, assemble_kmer, assemble_thread, True, True,
                      change_seed, kmer_limit, minimum_length_ratio, maximum_length_ratio, 1000, generate_scaffold),
                daemon=True)
            thread.start()
        if event == "-ASSEMBLE DONE-":
            thread = None
            print("INFO: All steps have been completed!")
    window.close()


if __name__ == "__main__":
    easy353_gui()
#     Brassicaceae
#  Apiaceae
