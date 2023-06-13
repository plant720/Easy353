#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:50 下午
# @Author     : zzhen
# @File       : easy353-gui.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.
import multiprocessing
import os
import platform
import re
import gzip
import shutil
import time
import psutil
import threading
from collections import defaultdict
from urllib.error import HTTPError
from urllib.request import urlretrieve
import PySimpleGUI as sg
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor

from Easy353Lib.assemble import assemble_reads_to_seq
# import my python file
from Easy353Lib.filter import get_file_lst, dnaseq2int, filter_read_for_ref
from build_database import generate_download_info, detect_classification
from build_database import generate_fasta_path, generate_gene_file


def make_ref_kmer_dict(ref_kmer_dict: defaultdict, ref_path: str, k_size: int):
    # all the ref files
    # 对ref_path_lst进行sort 保证顺序
    ref_path_lst = sorted(get_file_lst(ref_path))
    seq_num_count = 0
    # 设置kmer_mask
    kmer_mask = (1 << (k_size << 1)) - 1
    # 获取参考文件列表长度
    ref_path_lst_len = len(ref_path_lst)
    for ref_index in range(0, ref_path_lst_len):
        ref_file = open(ref_path_lst[ref_index], 'rt')
        # 取出文件第一行
        next(ref_file)
        for line in ref_file:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = next(ref_file, None)
            # Preserve the numbers and letters in the seq sequence
            # Remove line breaks
            ref_seq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            # 基因数量+1
            seq_num_count += 1
            ref_bin_seq, seq_len = dnaseq2int(ref_seq)
            # 对ref序列进行反向互补
            ref_bin_seqs_lst = [ref_bin_seq, dnaseq2int(ref_seq, True)]
            for rc_flag, ref_bin_seq in enumerate(ref_bin_seqs_lst):
                bin_len = 22 + ref_path_lst_len
                for seq_index in range(seq_len - k_size + 1):
                    # 从前向后算
                    # 正反链信息 occurrence 文件信息
                    # 第1位：辅助位 第2位：正反链 第3-22位：occurrence 第23-(22+len(ref_path_lst)): 文件信息
                    bin_kmer_info = (1 << (bin_len - 1)) + (rc_flag << (bin_len - 2))
                    # 从ref_bin_seq中获取指定长度的kmer
                    bin_kmer = (ref_bin_seq >> (seq_index << 1)) & kmer_mask
                    if bin_kmer in ref_kmer_dict:
                        bin_kmer_info = ref_kmer_dict[bin_kmer]
                    # occurrence +1
                    bin_kmer_info += (1 << ref_path_lst_len)
                    # 设置文件位
                    bin_kmer_info |= (1 << (ref_path_lst_len - ref_index - 1))
                    ref_kmer_dict[bin_kmer] = bin_kmer_info
        print('INFO: Num. of Hashed Ref Seq: {}'.format(seq_num_count))


def reads_bin_filter(ref_path: str, ref_kmer_dict: defaultdict, k_size: int, s_size: int, fq_1: str, fq_2: str,
                     out_dir: str, p_id: int, p_count: int, filter_pair_read: bool = True):
    t_start, t_end, reads_count, ref_path_lst, file_dict = time.perf_counter(), 0, 0, sorted(get_file_lst(ref_path)), {}
    # detect whether the input is paired-end reads
    paired_reads = not fq_1 == fq_2
    # 设置输出文件
    for file_index, one_ref_path in enumerate(ref_path_lst):
        ref_file_name = os.path.splitext(os.path.basename(one_ref_path))[0]
        # out_fmt=1 输出是fastq out_fmt=0 输出是fasta
        # 将输出文件的句柄存放在dict中
        if not paired_reads:
            file_dict[file_index] = open(os.path.join(out_dir, ref_file_name + '.' + str(p_id + 1) + '.fasta'), 'w')
        else:
            # paired-end reads: the reads are written to two files
            file_dict[file_index] = open(os.path.join(out_dir, ref_file_name + '_R1' + '.' + str(p_id + 1) + '.fasta'),
                                         'w'), open(
                os.path.join(out_dir, ref_file_name + '_R2' + '.' + str(p_id + 1) + '.fasta'), 'w')
    # 判断文件是否是gz文件
    gz_file = True if fq_1.endswith('.gz') else False
    infile_1 = gzip.open(fq_1, 'rt') if gz_file else open(fq_1, 'rt')
    if paired_reads: infile_2 = gzip.open(fq_2, 'rt') if gz_file else open(fq_2, 'rt')
    for _ in infile_1:
        reads_count += 1
        tmp_rec1 = [_, next(infile_1, None), next(infile_1, None), next(infile_1, None)]
        if paired_reads:  tmp_rec2 = [next(infile_2, None), next(infile_2, None), next(infile_2, None),
                                      next(infile_2, None)]
        if reads_count % p_count == p_id:
            # write kmers into the files
            read_seqs_lst = [tmp_rec1[1], tmp_rec2[1]] if paired_reads and filter_pair_read else [tmp_rec1[1]]
            for file_index in filter_read_for_ref(ref_kmer_dict, k_size, s_size, read_seqs_lst, len(ref_path_lst),
                                                  filter_pair_read):
                if paired_reads:
                    file_dict[file_index][0].writelines(['>', tmp_rec1[0], tmp_rec1[1]])
                    file_dict[file_index][1].writelines(['>', tmp_rec2[0], tmp_rec2[1]])
                else:
                    file_dict[file_index].writelines(['>', tmp_rec1[0], tmp_rec1[1]])
            if (reads_count << int(paired_reads)) % 1000000 == 0:
                t_end = time.perf_counter()
                # 交换时间
                t_start, t_end = t_end, t_end - t_start
                print('INFO: Handled {} M reads. {:.2f} s/M reads'.format((reads_count << int(paired_reads)) // 1000000,
                                                                          t_end))
    # 关闭文件句柄
    for _ in file_dict.keys():
        if paired_reads:
            file_dict[_][0].close()
            file_dict[_][1].close()
        else:
            file_dict[_].close()
    infile_1.close()
    if paired_reads: infile_2.close()


# initialize the window
def window_init():
    icon = None
    sg.theme('Default 1')
    # Parameters
    # basic_para_frame 
    basic_para_frame = [
        # paired-end data
        [sg.Text('Paired fq file 1', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key='-1-', size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Input file with forward paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False),
         sg.FileBrowse("View", font=("Consolas", 13), target='-1-')],
        # paired-end data
        [sg.Text('Paired fq file 2', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key='-2-', size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Input file with reverse paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False),
         sg.FileBrowse("View", font=("Consolas", 13), target='-2-')],
        # classification
        [sg.Text('Taxonomy', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key="-classification-", size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Input the interested classification.", readonly=False)],
        # reference
        [sg.Text('Reference dir', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key="-reference-", size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Input the reference directory that stores the ref files.", readonly=False),
         sg.FolderBrowse("View", font=("Consolas", 13), target='-reference-')],
        # output
        [sg.Text('Output dir', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key='-o-', size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip='Output directory.',
                  readonly=False),
         sg.FolderBrowse("View", font=("Consolas", 13), target='-o-')],
    ]
    # general_para_frame
    general_para_frame = [
        [sg.Text('Exclude file', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key='-exclude_file-', size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip='The file documenting species that need to be excluded',
                  readonly=False),
         sg.FileBrowse("View", font=("Consolas", 13), target='-exclude_file-')],

        [sg.Text("Exclude", size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(key="-exclude-", size=(11, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Input species that need to be excluded", readonly=False)],

        # kmer for filter and assembly
        [sg.Text('Filtering Kmer', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(21, key='-fk-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Kmer setting for filtering reads. [Default:21]"),
         sg.Text('Assembly Kmer', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(31, key='-ak-', size=(8, 1), font=("Consolas", 13), expand_x=True,
                  tooltip="Kmer setting for assembling reads. [Default:31]"
                  )],

        # threads for filter and assembly
        [sg.Text('Filtering thread', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(1, key='-ft-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Threads setting for filtering reads. [Default:4]"),
         sg.Text('Assembly thread', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(10, key='-at-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Threads setting for assembling reads. [Default:10]"
                  )],
    ]
    advanced_para_frame = [
        # sliding_window
        [sg.Text('Step length', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(1, key='-s-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Step length of the sliding window on the reads. [Default:1]"
                  ),
         sg.Text('Change seed', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(32, key='-change_seed-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Times of changing seed. [Default:32]")
         ],
        #  kmer_limit
        [sg.Text('Kmer limit', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(4, key='-kmer_limit-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="Limit of kmer count. [Default:4]"
                  ),
         sg.Text('Min length ratio', size=(11, 1), justification='right', font=("Consolas", 13)),
         sg.Input(0.8, key='-minimum_length_ratio-', size=(8, 1), font=("Consolas", 13), expand_x=True,

                  tooltip="The minimum ratio of contig length to reference average length. [Default:0.8]")],
        # fpr and gdk
        [sg.Text('Filter pair read', size=(13, 1), justification='center', font=("Consolas", 13)),
         sg.Checkbox(text="", default=False,
                     key='-filter_pair_read-', size=(8, 1), font=("Consolas", 13), ),
         sg.Text('Get dynamic kmer', size=(13, 1), justification='center', font=("Consolas", 13)),
         sg.Checkbox(text="", default=False, key='-get_dynamic_kmer-', size=(8, 1), font=("Consolas", 13), )],
    ]

    # left_col
    left_col = [
        # parameter frame
        [sg.Frame('Basic', layout=basic_para_frame, expand_x=True, font=("Consolas", 13, 'bold'))],
        [sg.Frame('General', layout=general_para_frame, expand_x=True, font=("Consolas", 13, 'bold'))],
        [sg.Frame('Advanced', layout=advanced_para_frame, expand_x=True, font=("Consolas", 13, 'bold'))],
    ]
    # right_col
    right_col = [
        [sg.Multiline(size=(60, 20), font=("Consolas", 11), expand_x=True, expand_y=True,
                      key='-OUTPUT-', autoscroll=True, auto_refresh=True, reroute_stdout=True, reroute_stderr=True)],
        # button
        [sg.Button('Run', font=("Consolas", 13), auto_size_button=True),
         sg.Button('Reset', font=("Consolas", 13), auto_size_button=True),
         sg.Button('Close', font=("Consolas", 13), auto_size_button=True)],
    ]
    info_frame = [
        [sg.Text('Authors: zzhen', font=("Consolas", 12), justification='left'),
         sg.Text('Email: zzhen0302@163.com', font=("Consolas", 12), justification='left')],
    ]
    layout = [
        [sg.Push(), sg.Column(info_frame, element_justification='left'), sg.Push()],
        [sg.Pane([sg.Column(left_col, element_justification='c', expand_x=True, expand_y=True),
                  sg.Column(right_col, element_justification='c', expand_x=True, expand_y=True)],
                 orientation='horizontal', relief=sg.RELIEF_SUNKEN, key='-PANE-')],
    ]

    window = sg.Window(title='Easy353 v2.0.2', layout=layout, icon=icon, finalize=True)
    return window


def net_connected():
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        if s.connect_ex(("www.baidu.com", 443)) == 0 or s.connect_ex(("www.google.com", 443)) == 0:
            s.close()
            return True
    except ConnectionError:
        return False


# function used to get a fasta file based on an url
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
        except Exception as e:
            info = "ERROR: {}".format(e)
    return info


# download files based on the taxon
def download_species_thread(output_dir: str, classification: list, window: sg.Window, max_threads: int = 8):
    print("INFO: Start of generating reference...")
    if output_dir is not None:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    print("INFO: Get species information that needs to be downloaded")
    # Select the species to be downloaded based on the taxon
    down_spec_info = generate_download_info(detect_classification(classification))
    # Detect existing file names
    existed_files = [os.path.basename(i) for i in generate_fasta_path(output_dir)]
    # Get files that have not been downloaded yet
    down_spec_info = [i for i in down_spec_info if i["Fasta file name"] not in existed_files]
    print("INFO: Download species data")
    print("INFO: Total {} species need to be downloaded".format(len(down_spec_info)))
    count = 1
    # use concurrent.features.ThreadPoolExecutor to speed up the download
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        tasks = []
        for spec_info in down_spec_info:
            tasks.append(executor.submit(download_files, spec_info, output_dir))
        count = 0
        for future in as_completed(tasks):
            count += 1
            if future.result() is not None:
                print(future.result())
            if count % 10 == 0:
                print("INFO: {} / {} has been downloaded".format(count, len(down_spec_info)))
            if count == len(down_spec_info):
                print("INFO: All species have been downloaded")
    print("INFO: Download finished!")
    # communicating to the main process
    window.write_event_value('-DOWNLOAD DONE-', "-DOWNLOAD DONE-")


# Generate reference sequences based on downloaded 353 data and species to be excluded
def generate_ref_thread(input_dir: str, output_dir: str, exclude: list, exclude_file: str, window: sg.Window):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
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
    file_path_list = generate_fasta_path(input_dir)
    generate_gene_file(_file_path_list_=file_path_list,
                       _output_dir_=output_dir, _exclude_species_=exclude_species)
    print("INFO: End of generating reference!")
    window.write_event_value("-GENERATE DONE-", "-GENERATE DONE-")


def easy353_thread(window: sg.Window, fq_file_1: str, fq_file_2: str, reference: str, out_dir: str,
                   filter_kmer: int = 21, filter_thread: int = 1, assemble_kmer=31, assemble_thread: int = 10,
                   step_length: int = 1, change_seed: int = 32, kmer_limit: int = 4, filter_pair_read: bool = False,
                   get_dynamic_kmer: bool = False):
    filter_out = os.path.join(out_dir, 'filter_out')
    if not os.path.isdir(filter_out):
        os.makedirs(filter_out)
    assemble_out = os.path.join(out_dir, 'assemble_out')
    if not os.path.isdir(assemble_out):
        os.makedirs(assemble_out)
    print("INFO: Easy353 is running ")
    print("INFO: 1. Build reference hash table ")

    t_hash_start = time.perf_counter()
    ref_kmer_dict = defaultdict(int)
    make_ref_kmer_dict(ref_kmer_dict, reference, filter_kmer)
    print('INFO: Hash dictionary has been made.')
    t_hash_end = time.perf_counter()
    print('INFO: Time used for building hash table: {:.2f} s'.format(t_hash_end - t_hash_start))
    print('INFO: The memory usage of the current process is: {:.2f} GB'.
          format(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024))
    print("INFO: 2. Read filtering ")
    print('INFO: Filter reads from fq_files based on hash table with {} thread'.format(filter_thread))
    # 处理输入的测序文件
    if not fq_file_1 and not fq_file_2:
        raise RuntimeError('There were no sequencing files inputted.')

    if filter_thread == 1:
        # when filter_thread is 1
        reads_bin_filter(reference, ref_kmer_dict, filter_kmer, step_length,
                         fq_file_1, fq_file_2,
                         filter_out, 0, filter_thread, filter_pair_read)
    else:
        process_list = []
        for p_id in range(filter_thread):
            p = multiprocessing.Process(target=reads_bin_filter,
                                        args=(reference, ref_kmer_dict, filter_kmer, step_length,
                                              fq_file_1, fq_file_2,
                                              filter_out, p_id, filter_thread, filter_pair_read,))
            process_list.append(p)
        for p in process_list: p.start()
        for p in process_list: p.join()

    # 合并多进程产生的文件
    ref_path_lst = get_file_lst(reference)
    paired_reads = True if fq_file_1 != fq_file_2 else False
    for ref_file in ref_path_lst:
        gene_name = os.path.splitext(os.path.basename(ref_file))[0]
        prefixes = [gene_name + "_R1", gene_name + "_R2"] if paired_reads else [gene_name]
        for prefix in prefixes:
            with open(os.path.join(filter_out, prefix + ".fasta"), 'wt') as merge_file:
                for file_name in [os.path.join(filter_out, prefix + "." + str(i + 1) + ".fasta") for i in
                                  range(filter_thread)]:
                    with open(file_name, 'rt') as infile:
                        shutil.copyfileobj(infile, merge_file)
                    os.remove(file_name)

    t_filter_end = time.perf_counter()
    print('INFO: Time used for filter: {:.2f} s'.format(t_filter_end - t_hash_end))
    print('INFO: Total time used for filter: {:.2f} s'.format(t_filter_end - t_hash_start))

    t_assemble_start = time.perf_counter()
    print('INFO: 3. Read Assembly ')
    print('INFO: {} genes will be assembled with using {} threads.'.format(len(ref_path_lst), assemble_thread))
    if assemble_thread == 1:
        # when filter_thread is 1
        for ref_file_name in ref_path_lst:
            gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
            fa_1 = os.path.join(filter_out, gene_name + "_R1.fasta") if paired_reads else os.path.join(
                filter_out, gene_name + ".fasta")
            fa_2 = os.path.join(filter_out, gene_name + "_R2.fasta") if paired_reads else fa_1
            assemble_reads_to_seq(gene_name, fa_1, fa_2, ref_file_name, assemble_kmer, assemble_out,
                                  kmer_limit, change_seed, get_dynamic_kmer)
    else:
        with ProcessPoolExecutor(max_workers=min(assemble_thread, len(ref_path_lst))) as executor:
            tasks = []
            for ref_file_name in ref_path_lst:
                gene_name = os.path.splitext(os.path.basename(ref_file_name))[0]
                fa_1 = os.path.join(filter_out, gene_name + "_R1.fasta") if paired_reads else os.path.join(
                    filter_out, gene_name + ".fasta")
                fa_2 = os.path.join(filter_out, gene_name + "_R2.fasta") if paired_reads else fa_1
                tasks.append(executor.submit(assemble_reads_to_seq, gene_name, fa_1, fa_2, ref_file_name,
                                             assemble_kmer, assemble_out, kmer_limit, change_seed, get_dynamic_kmer))
            count = 0
            for future in as_completed(tasks):
                count += 1
                if count % 10 == 0:
                    print("INFO: {} / {} genes have been assembled!".format(count, len(ref_path_lst)))
                if count == len(ref_path_lst):
                    print("INFO: All genes have been assembled!".format(count, len(ref_path_lst)))
    t_assemble_end = time.perf_counter()
    print('INFO: Total time used for assembly: {:.2f} s'.format(t_assemble_end - t_assemble_start))
    print('INFO: All jobs have been done!')
    print('INFO: Total time used: {:.2f} s'.format(t_assemble_end - t_hash_start))
    window.write_event_value("-EASY353 DONE-", "-EASY353 DONE-")


def easy353_gui():
    window = window_init()
    fq_file_1, fq_file_2, output_dir = None, None, None
    filter_kmer, assemble_kmer, filter_thread, assemble_thread = 21, 31, 1, 10
    step_length, kmer_limit, minimum_length_ratio, change_seed = 1, 4, 0.8, 32
    filter_pair_read, get_dynamic_kmer = False, False

    # get_ref parameter
    reference, classification, exclude_file, exclude = None, None, None, None
    # Used to prevent double click events
    thread = None

    # set the output_dir
    species_dir, ref_dir, easy353_out = None, None, None

    while True:
        # filter and assemble parameter
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Close':  # if user closes window or clicks cancel
            break
        if event == 'Reset' and not thread:
            window['-1-'].update('')
            window['-2-'].update('')
            window['-o-'].update('')
            window['-classification-'].update('')
            window['-reference-'].update('')

            window['-exclude_file-'].update('')
            window['-exclude-'].update('')
            window['-fk-'].update(21)
            window['-ak-'].update(31)
            window['-ft-'].update(1)
            window['-at-'].update(10)

            window['-s-'].update(1)
            window['-change_seed-'].update(32)
            window['-kmer_limit-'].update(2)
            window['-minimum_length_ratio-'].update(0.8)
            window['-filter_pair_read-'].update(False)
            window['-get_dynamic_kmer-'].update(False)

            fq_file_1, fq_file_2, output_dir = None, None, None
            filter_kmer, assemble_kmer, filter_thread, assemble_thread = 21, 31, 1, 10
            step_length, kmer_limit, minimum_length_ratio, change_seed = 1, 4, 0.8, 32
            filter_pair_read, get_dynamic_kmer = False, False

            # get_ref parameter
            reference, classification, exclude_file, exclude = None, None, None, None
            # output dir
            species_dir, ref_dir, easy353_out = None, None, None

        if event == 'Run' and not thread:
            # basic parameters
            if values["-1-"]:
                fq_file_1 = values["-1-"]
            if values["-2-"]:
                fq_file_2 = values["-2-"]
            if not fq_file_1 and not fq_file_2:
                sg.Popup("Please input fastq file(s)!",
                         title='Info', keep_on_top=True, font=("Consolas", 13))
                continue
            if fq_file_1 and not fq_file_2:
                fq_file_2 = fq_file_1
            if fq_file_2 and not fq_file_1:
                fq_file_1 = fq_file_2
            if values["-o-"]:
                output_dir = values["-o-"]
            else:
                sg.Popup("Please input an output directory!",
                         title='Info', keep_on_top=True, font=("Consolas", 13))
                continue
            if values["-classification-"]:
                # classification
                classification = values["-classification-"]
                classification = re.split(r"[,，、 ]", classification)
            if values["-reference-"]:
                reference = values["-reference-"]
            if not values['-reference-'] and not values['-classification-']:
                sg.Popup("Please input a reference directory path or a classification!", title='Info', keep_on_top=True,
                         font=("Consolas", 13))
                continue

            if values["-exclude-"]:
                exclude = values["-exclude-"]
                exclude = re.split(r"[,，、 ]", exclude)
            if values["-exclude_file-"]:
                exclude_file = values["-exclude_file-"]
            if values["-fk-"]:
                try:
                    filter_kmer = int(values["-fk-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of filtering.kmer should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-ak-"]:
                try:
                    assemble_kmer = int(values["-ak-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of assembly.kmer should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-ft-"]:
                try:
                    filter_thread = int(values["-ft-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of filtering.thread should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-at-"]:
                try:
                    assemble_thread = int(values["-at-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of assembly.thread should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-s-"]:
                try:
                    step_length = int(values["-s-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of step.length should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-change_seed-"]:
                try:
                    change_seed = int(values["-change_seed-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of change.seed should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-kmer_limit-"]:
                try:
                    kmer_limit = int(values["-kmer_limit-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of kmer.limit should be integer!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-minimum_length_ratio-"]:
                try:
                    minimum_length_ratio = float(values["-minimum_length_ratio-"])
                except (ValueError, TypeError):
                    sg.Popup("The value of min.length.ratio should be float!",
                             title='Info', keep_on_top=True, font=("Consolas", 13))
                    continue
            if values["-filter_pair_read-"]:
                filter_pair_read = True
            if values['-get_dynamic_kmer-']:
                get_dynamic_kmer = True

            window.write_event_value("-PARAMETERS DONE-", "-PARAMETERS DONE-")
        if event == '-PARAMETERS DONE-':
            if reference:
                print("INFO: The files under the reference directory will be used as reference!")
                if classification:
                    print("INFO: The classification will be ignored!")
                window.write_event_value("-GENERATE DONE-", "-GENERATE DONE-")
            if classification and not reference:
                print(
                    "INFO: The 353species files will be downloaded from treeoflife.kew.org based on the specific-taxonomy, which will be used as reference!")
                # download the AGS sequences from Kew Tree of Life Explorer
                species_dir = os.path.join(output_dir, "353species_files")
                # check the net connect
                if not net_connected():
                    sg.Popup("The network is not connected,please check it! The software will quit later!",
                             title='Warning', keep_on_top=True, font=("Consolas", 13))
                    window.close()
                thread = threading.Thread(
                    target=download_species_thread,
                    args=(species_dir, classification, window),
                    daemon=True)
                thread.start()
        if event == '-DOWNLOAD DONE-':
            thread.join(timeout=0)
            # After downloading the AGS sequences, generate the reference files
            # The variable thread has changed
            ref_dir = os.path.join(output_dir, "353genes_ref")
            # set the reference as ref_dir
            reference = ref_dir
            # if None file under species_dir
            if len(os.listdir(species_dir)) == 0:
                sg.Popup("No file was downloaded, please check whether the taxonomy: " + ",".join(classification)
                         + " is right! The software will quit later!", title='Error',
                         keep_on_top=True, font=("Consolas", 13))
                window.close()
            thread = threading.Thread(
                target=generate_ref_thread,
                args=(species_dir, ref_dir, exclude, exclude_file, window),
                daemon=True)
            thread.start()
        if event == '-GENERATE DONE-':
            if thread:
                thread.join(timeout=0)
            # After building the reference sequences, enter the new thread to do read filtering
            easy353_out = os.path.join(output_dir, "easy353_out")
            thread = threading.Thread(
                target=easy353_thread,
                args=(window, fq_file_1, fq_file_2, reference, easy353_out, filter_kmer, filter_thread, assemble_kmer,
                      assemble_thread, step_length, change_seed, kmer_limit, filter_pair_read, get_dynamic_kmer),
                daemon=True)
            thread.start()
        if event == "-ASSEMBLE DONE-":
            thread.join(timeout=0)
            # set thread None
            thread = None
            print("INFO: All steps have been completed!")

    window.close()


if __name__ == "__main__":
    if platform.system() in ("Linux", "Darwin"):
        multiprocessing.set_start_method('fork')
    else:
        multiprocessing.set_start_method('spawn')
    easy353_gui()
