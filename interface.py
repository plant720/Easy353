#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:50 下午
# @Author     : zzhen
# @File       : interface.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import PySimpleGUI as sg
from src.assemble import assemble_flow
from src.filter import filter_flow
from src.build_database import *
import re


def main_window():
    sg.theme('Default 1')
    # Parameters
    # basic_para_frame
    basic_para_frame = [
        # sequence_data
        [sg.Text('Unpaired_fq_file', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-q-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with unpaired (single-end) reads.",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-q-')],

        # paired-end data
        [sg.Text('Fq_file_1', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with forward paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-1-')],
        # paired-end data
        [sg.Text('Fq_file_2', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input file with reverse paired-end reads (*.fq/.gz/.tar.gz).",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Times", 12), target='-2-')],

        # reference
        [sg.Text('Reference', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-r-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip='Input a directory with references.',
                  readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Times", 12), target='-r-')],

        # reads directory
        [sg.Text('Reads_dir', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-read-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input directory with filtered reads if you want to assemble from reads.",
                  readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Times", 12), target='-read-')],

        # reads directory
        [sg.Text('Species_dir', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-species_dir-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input directory with species",
                  readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Times", 12), target='-species_dir-')],

        # output
        [sg.Text('Output_dir', size=(14, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-o-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip='Output directory.',
                  readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Times", 12), target='-o-')],
    ]
    # general_para_frame
    general_para_frame = [
        # classification
        [sg.Text('Classification', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(key="-classification-", size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input the interested classification.", readonly=False, background_color='#BFE9FF')],

        [sg.Text('Exclude_file', size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(key='-exclude_file-', size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip='The file documenting species that need to be excluded',
                  readonly=False, background_color='#BFE9FF'),
         sg.FileBrowse("File", font=("Times", 12), target='-exclude_file-')],

        [sg.Text("Exclude", size=(18, 1), justification='center', font=("Times", 12)),
         sg.Input(key="-exclude-", size=(12, 1), font=("Times", 12), expand_x=True,
                  tooltip="Input species that need to be excluded", readonly=False,
                  background_color='#BFE9FF')],

        # kmer for filter and assembly
        [sg.Text('K-mer for filter', size=(18, 1), justification='right', font=("Times", 12)),
         sg.Input(31, key='-k1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Kmer setting for filtering reads. [Default:31]"),
         sg.Text('K-mer for assembly', size=(18, 1), justification='right', font=("Times", 12)),
         sg.Input(41, key='-k2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Kmer setting for assembling reads. [Default:41]"
                  )],

        # threads for filter and assembly
        [sg.Text('Threads for filtering', size=(18, 1), justification='right', font=("Times", 12)),
         sg.Input(1, key='-t1-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Threads setting for filtering reads. [Default:4]"),
         sg.Text('Threads for assembly', size=(18, 1), justification='right', font=("Times", 12)),
         sg.Input(1, key='-t2-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Threads setting for assembling reads. [Default:4]"
                  )],
    ]
    advanced_para_frame = [
        # sliding_window
        [sg.Text('Step_length', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(1, key='-s-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Step length of the sliding window on the reads. [Default:1]"
                  ),
         sg.Text('Ref_ratio', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(1.0, key='-reference_ratio-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The percentage of the reference sequence used. [Default:1.0]")],
        # change_seed and limit_count
        [sg.Text('Change_seed', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(32, key='-change_seed-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Times of changing seed. [Default:32]"),

         sg.Text('Kmer_limit', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(2, key='-kmer_limit-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="Limit of kmer count. [Default:2]"
                  )],
        # min_percent_length and max_percent_length
        [sg.Text('Min_length_ratio', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(1.0, key='-minimum_length_ratio-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The minimum ratio of contig length to reference average length. [Default:1.0]"),

         sg.Text('Max_length_ratio', size=(14, 1), justification='right', font=("Times", 12)),
         sg.Input(2.0, key='-maximum_length_ratio-', size=(12, 1), font=("Times", 12), expand_x=True,
                  background_color='#BFE9FF',
                  tooltip="The maximum ratio of contig length to reference average length. [Default:2.0]")
         ],

        # reverse_complement
        [sg.Text('Ref_reverse_complement', size=(21, 1), justification='right', font=("Times", 12)),
         sg.Checkbox(text="True", default=True,
                     key='-reverse_complement-', size=(12, 1), font=("Times", 12), ),

         sg.Text('Generate_scaffold', size=(15, 1), justification='right', font=("Times", 12)),
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
        [sg.Text("Function", size=(12, 1), justification='center', font=("Times", 12, 'bold')),
         sg.Checkbox("get_ref", default=True, key='-reference-', font=("Times", 12)),
         sg.Checkbox("filter", default=False, key='-filter-', font=("Times", 12)),
         sg.Checkbox("assembly", default=False, key='-assembly-', font=("Times", 12))
         ],
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


def easy353_gui():
    window = main_window()
    running_flag = False
    while True:
        # filter and assemble parameter
        fastq_files, _paired_reads_ = tuple(), False
        fq_file_1, fq_file_2, unpaired_fq_file, reads_input_dir, reference, output_dir \
            = None, None, None, None, None, None
        filter_kmer, assemble_kmer, filter_thread, assemble_thread, function_mode = 31, 41, 1, 1, 0
        step_length, kmer_limit, minimum_length_ratio, maximum_length_ratio, change_seed, reference_ratio \
            = 1, 2, 1.0, 2.0, 32, 1.0
        # fast 即 ref_reverse_complement
        generate_scaffold, fast = False, False
        # get_ref parameter
        species_dir, classification, exclude_file, exclude = None, None, None, None
        ref_flag, filter_flag, assembly_flag = False, False, False

        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Close':  # if user closes window or clicks cancel
            break
        if event == 'Reset' and not running_flag:
            window['-1-'].update('')
            window['-2-'].update('')
            window['-q-'].update('')
            window['-read-'].update('')
            window['-r-'].update('')
            window['-species_dir-'].update('')
            window['-o-'].update('')

            window['-classification-'].update('')
            window['-exclude_file-'].update('')
            window['-exclude-'].update('')
            window['-k1-'].update(31)
            window['-k2-'].update(41)
            window['-t1-'].update(1)
            window['-t2-'].update(1)
            window['-s-'].update(1)
            window['-reference_ratio-'].update(1.0)

            window['-change_seed-'].update(32)
            window['-kmer_limit-'].update(2)
            window['-minimum_length_ratio-'].update(1.0)
            window['-maximum_length_ratio-'].update(2.0)

            window['-generate_scaffold-'].update(False)
            window['-reverse_complement-'].update(True)

            window['-reference-'].update(True)
            window['-filter-'].update(False)
            window['-assembly-'].update(False)

        if event == 'Run' and not running_flag:
            # 设定running_flag为True，防止重复点击Run
            running_flag = True
            if values["-reference-"]:
                ref_flag = True
            if values["-filter-"]:
                filter_flag = True
            if values["-assembly-"]:
                assembly_flag = True

            if ref_flag:
                if values["-species_dir-"]:
                    species_dir = values["-species_dir-"]
                if values["-classification-"]:
                    classification = values["-classification-"]
                    classification = re.split(r"[,，、 ]", classification)
                    print(classification)
                if not classification and not species_dir:
                    sg.Popup("Please input a directory or a classification!", title='Info', keep_on_top=True,
                             font=("Times", 12))
                if values["-exclude-"]:
                    exclude = values["-exclude-"]
                    exclude = re.split(r"[,，、 ]", exclude)
                if values["-exclude_file-"]:
                    exclude_file = values["-exclude_file-"]

            if filter_flag:
                # 基本参数部分
                if values["-1-"]:
                    fq_file_1 = values["-1-"]
                if values["-2-"]:
                    fq_file_2 = values["-2-"]
                if values["-q-"]:
                    unpaired_fq_file = values["-q-"]
                if unpaired_fq_file or (fq_file_1 and fq_file_2):
                    if unpaired_fq_file:
                        fastq_files = ([unpaired_fq_file],)
                    if fq_file_1 and fq_file_2:
                        fastq_files = ([fq_file_1], [fq_file_2])
                        _paired_reads_ = True
                else:
                    sg.Popup("Please input fastq file(s)!",
                             title='Info', keep_on_top=True, font=("Times", 12))
                    continue
            if values["-r-"]:
                reference = values["-r-"]
            if values["-o-"]:
                output_dir = values["-o-"]
            if assembly_flag:
                if not filter_flag:
                    if values["-read-"]:
                        reads_input_dir = values["-read-"]
                    else:
                        sg.Popup("Please input reads directory!",
                                 title='Info', keep_on_top=True, font=("Times", 12))
                        continue
                else:
                    reads_input_dir = output_dir

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
            if values["-reference_ratio-"]:
                reference_ratio = float(values["-reference_ratio-"])
                if reference_ratio > 1.0 or reference_ratio < 0:
                    sg.Popup("Reference ratio should be between 0 and 1!",
                             title='Info', keep_on_top=True, font=("Times", 12))
                    continue
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
            # fq_file_1 = "/Users/zzhen/Desktop/ara.1.fq"
            # fq_file_2 = "/Users/zzhen/Desktop/ara.2.fq"
            # reference = "/Users/zzhen/Desktop/new"
            # output_dir = "/Users/zzhen/Desktop/test_a"
            # filter_thread = 1
            # assembly_thread = 1
            # fastq_files = ([fq_file_1], [fq_file_2])
            # _paired_reads_ = True

            if ref_flag:
                print("INFO: Start of generating reference...")
                # 处理输入和输出文件夹
                exclude_species = []
                if species_dir is None and output_dir is None:
                    print("Input and output directory cannot both be empty!")
                    exit(-1)
                if output_dir is not None:
                    if not os.path.isdir(output_dir):
                        os.makedirs(output_dir)
                    if species_dir is None:
                        species_dir = output_dir
                if species_dir is not None and output_dir is None:
                    output_dir = species_dir

                # 当需要下载文件时
                if classification is not None:
                    print("INFO: Get species information that needs to be downloaded")
                    # 根据传入的分类信息需要下载的物种信息 是元素为dict的list
                    down_spec_info = generate_download_info(detect_classification(classification))
                    # 检测已经存在的文件名
                    existed_files = [i.split("/")[1] for i in generate_fasta_path(species_dir)]
                    # 获取还未下载的文件
                    down_spec_info = [i for i in down_spec_info if i not in existed_files]
                    # 对于多参数函数，如果我们只想对它的一个参数在多进程任务中依次取可迭代对象中各个值，其他参数固定，
                    # 可以使用偏函数构造出单参数函数
                    print("INFO: Download species data")
                    # if _thread_ > 1:
                    #     with Pool(processes=_thread_) as pool:
                    #         result = list(tqdm(pool.imap(partial(download_dict, input_dir=input_dir), down_spec_info),
                    #                            total=len(down_spec_info), desc="download files"))
                    #         print(result)
                    # else:
                    bar_format = '{desc}{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt}'
                    with tqdm(total=len(down_spec_info), desc="download files:", bar_format=bar_format) as _tqdm:
                        for _spec_info_ in down_spec_info:
                            download_dict(_spec_info_, species_dir)
                            _tqdm.update(1)
                if exclude is not None:
                    exclude_species.extend(exclude)
                if exclude_file is not None:
                    if not os.path.isfile(exclude_file):
                        print("The exclude file does not exist!")
                    else:
                        with open(exclude_file, "r") as _file_:
                            _exclude_ = _file_.readlines()
                            _exclude_ = [x.strip() for x in _exclude_ if x.strip()]
                            exclude_species.extend(_exclude_)
                exclude_species = [species.capitalize().replace(" ", "_") for species in exclude_species]

                print("INFO: Generating species data into reference files")
                # 将input_dir下的文件和 output_dir中的文件合并
                input_file_list = generate_fasta_path(species_dir)
                output_file_list = generate_fasta_path(output_dir)
                file_path_list = list((set(input_file_list).union(set(output_file_list))))
                generate_gene_file(_file_path_list_=file_path_list,
                                   _output_dir_=os.path.join(output_dir, "353gene"), _exclude_species_=exclude_species)
                print("INFO: End of generating reference!")

            if filter_flag:
                print("INFO: Start of reads filtering...")
                filter_flow(_read_data_tuple_=fastq_files, _out_dir_=output_dir,
                            _reference_path_=reference, _kmer_size_=filter_kmer,
                            _step_size_=step_length, _ref_reverse_complement_=fast,
                            _paired_reads_=_paired_reads_, _thread_for_filter_=filter_thread, _print_=False)
                print("INFO: End of reads filtering!")
                # 当点击run后不允许再次运行
            if assembly_flag:
                print("INFO: Start of reads assembly...")
                assemble_flow(_input_read_path_=reads_input_dir, _out_dir_=output_dir, _ref_path_=reference,
                              _assemble_kmer_size_=assemble_kmer, _assemble_thread_=assemble_thread,
                              _ref_reverse_complement_=True, _pos_=True,
                              _change_seed_=change_seed, _kmer_limit_count_=kmer_limit,
                              _min_percent_length_=minimum_length_ratio,
                              _max_percent_length_=maximum_length_ratio,
                              _iteration_=1000, _write_scaffold_=generate_scaffold)
                print("INFO: End of reads assembly!")
            # 当点击run后不允许再次运行
            running_flag = False
    window.close()


if __name__ == "__main__":
    easy353_gui()
