#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/4/6 2:50 下午
# @Author     : zzhen
# @File       : interface.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import PySimpleGUI as sgui


if __name__ == "__main__":

    sgui.theme('Default 1')  # Add a touch of color
    # All the stuff inside your window.

    # data_frame

    data_frame = [
        # reference
        [sgui.Text('参考序列:', size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Input(key='-r-', size=(20, 1), font=("宋体", 12), expand_x=False,
                    tooltip='选择一个参考序列文件',
                    readonly=False, background_color='#BAD9BC'),
         sgui.FolderBrowse("Folder", font=("Times", 12), target='-r-')],

        # sequence_data
        [sgui.Text('单端测序数据:', size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Input(key='-q-', size=(20, 1), font=("宋体", 12), expand_x=False,
                    tooltip="输入单端测试数据",
                    readonly=False, background_color='#BAD9BC'),
         sgui.FileBrowse("File", font=("Times", 12), target='-q-')],

        # paired-end data
        [sgui.Text('双端测序数据PE-1:', size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Input(key='-q1-', size=(20, 1), font=("宋体", 12), expand_x=False,
                    tooltip="输入双端测试数据-1",
                    readonly=False, background_color='#BAD9BC'),
         sgui.FileBrowse("File", font=("Times", 12), target='-q1-')],
        # paired-end data
        [sgui.Text('双端测序数据PE:', size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Input(key='-q2-', size=(20, 1), font=("宋体", 12), expand_x=False,
                    tooltip="输入双端测试数据-2",
                    readonly=False, background_color='#BAD9BC'),
         sgui.FileBrowse("File", font=("Times", 12), target='-q2-')],

        # output
        [sgui.Text('输出文件夹:', size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Input(key='-o-', size=(20, 1), font=("宋体", 12), expand_x=False,
                    tooltip='指定一个文件夹用于输出数据',
                    readonly=False, background_color='#BAD9BC'),
         sgui.FolderBrowse("Folder", font=("Times", 12), target='-o-')],

        [sgui.Text("程序功能选择", size=(18, 1), justification='center', font=("宋体", 12)),
         sgui.Radio('执行所有功能', "function", default=True, key="-all-"),
         sgui.Radio('根据参考序列过滤测序数据', "function", default=False, key="-filter-"),
         sgui.Radio('对过滤的测序数据重新过滤', "function", default=False, key="-refilter-"),
         ],
    ]

    # option_frame
    option_frame = [

        # kmer for filter and assembly
        [sgui.Text('K-mer for filter:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(29, key='-k1-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the k-mer for filter [default=29]"),
         sgui.Text('K-mer for assembly:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(31, key='-k2-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the k-mer for assembly [default=31]"
                    )],
        # threads for filter and assembly
        [sgui.Text('Threads for filter:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(4, key='-t1-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The number of threads for filter [default=4]"),
         sgui.Text('Threads for assembly:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(8, key='-t2-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The number of threads for assembly [default=8]"
                    )],
        # sliding_window
        [sgui.Text('Sliding window:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(1, key='-s-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The size of the sliding window [default=1]"
                    )],
    ]

    default_frame = [
        # change_seed and limit_count
        [sgui.Text('Change_seeds:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(32, key='-change_seed-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the change_seed [default=32]"),

         sgui.Text('Limit_count:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(2, key='-limit_count-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the limit_count [default=2]"
                    )],
        # min_percent_length and max_percent_length
        [sgui.Text('Min_percent_length:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(0.5, key='-min_percent_length-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the min_percent_length [default=0.5]"),

         sgui.Text('Max_percent_length:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(1.5, key='-max_percent_length-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the max_percent_length [default=1.5]")
         ],

        # reads_length
        [sgui.Text('Reads_length:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(150, key='-reads_length-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the reads_length [default=150]"),
         # contigs
         sgui.Text('Contigs:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Input(10, key='-contigs-', size=(20, 1), font=("Times", 12), expand_x=False,
                    background_color='#BFE9FF',
                    tooltip="The value of the contigs [default=100]"),
         ],
        # reverse_complete
        [sgui.Text('Reverse_complete:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Checkbox("参考序列反向互补", default=True, key='-reverse_complete-', size=(20, 1), font=("Times", 12), ),
         # python3
         sgui.Text('Python3:', size=(18, 1), justification='right', font=("Times", 12)),
         sgui.Checkbox("python3", default=True, key='-python3-', size=(20, 1), font=("Times", 12), ),
         ],
    ]

    info_frame = [
        [sgui.Text('EasyGene', justification="center", font=("Times", 24, "bold"))
         ],
        [sgui.Text('Authors: Zhang Zhen', font=("Times", 12)),
         sgui.Text('Email: zzhen0302@163.com', font=("Times", 12))],
        [sgui.Text('Version: 1.0', font=("Times", 12))]
    ]

    layout = [
        info_frame,
        data_frame, option_frame, default_frame,
        [sgui.Text('程序运行记录', justification='center', font=("宋体", 12))],
        [sgui.Output(size=(100, 20), font=("Arial", 10))],
        # 按钮
        [sgui.Button('Run', size=(8, 1), font=("Times", 12)),
         sgui.Button('Exit', size=(8, 1), font=("Times", 12))]
    ]

    window = sgui.Window('基因生成器', layout)

    while True:
        event, values = window.read()
        if event == sgui.WIN_CLOSED or event == 'Exit':  # if user closes window or clicks cancel
            break
        if event == 'Run':
            # 基本参数部分
            if values["-f-"]:
                word_file = values["-f-"]
            else:
                sgui.Popup("必须传入参考序列文件",
                           title='Info', keep_on_top=True, font=("宋体", 12))
                continue
            if values["-2-"]:
                excel_file = values["-2-"]
            else:
                sgui.Popup("必须传入一个excel文件",
                           title='Info', keep_on_top=True, font=("宋体", 12))
                continue
            if values["-3-"]:
                output_dir = values["-3-"]
            else:
                sgui.Popup("必须传入一个输出文件夹",
                           title='Info', keep_on_top=True, font=("宋体", 12))
                continue

    window.close()