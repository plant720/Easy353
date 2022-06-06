from pickle import TRUE
import assemble
from collections import defaultdict
import time
import os, shutil
import gzip
from Bio import SeqIO
import sys
from Bio.Align.Applications import MuscleCommandline
import argparse
import gc
from numpy import average
import subprocess, threading
from multiprocessing import Process
import multiprocessing

pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=''' Miner YY ''')
pars.add_argument('-f', metavar='<str>', type=str, help='''input fasta files.''', required=False,
                  default='./refs/Apiale/5922.fasta')
pars.add_argument('-q', metavar='<str>', type=str, help='''input fastq files.''', required=False, nargs="+")
pars.add_argument('-q1', metavar='<str>', type=str, help='''input fastq -1 files.''', required=False, nargs="+")
pars.add_argument('-q2', metavar='<str>', type=str, help='''input fastq -2 files.''', required=False, nargs="+")
pars.add_argument('-k1', metavar='<int>', type=int, help='''kmer of filter''', default=31)
pars.add_argument('-k2', metavar='<int>', type=int, help='''kmer of assemble''', default=41)
pars.add_argument('-s', metavar='<int>', type=int, help='''the length of the sliding window on the reads''', default=1)
pars.add_argument("-t0", metavar='<int>', type=int, help="Tthread of read file", default=4)
pars.add_argument("-t", metavar='<int>', type=int, help="Tthread", default=8)
pars.add_argument('-o', metavar='<str>', type=str, help='''out dir.''', required=False, default='./results_AP_353_3/')
pars.add_argument('-mode', metavar='<int>', type=int, help='''0:all, 1:filter, 2: re-filter, 3: assemble ''',
                  required=False, default=3)
pars.add_argument('-change_seed', metavar='<int>', type=int, help='''times of change seed''', required=False,
                  default=32)
pars.add_argument('-limit_count', metavar='<int>', type=int, help='''limit of kmer count''', required=False, default=-1)
pars.add_argument('-limit_min_length', metavar='<float>', type=float, help='''limit of contig length''', required=False,
                  default=0.5)
pars.add_argument('-limit_max_length', metavar='<float>', type=float, help='''limit of contig length''', required=False,
                  default=2)
pars.add_argument('-list', metavar='<str>', type=str, help='''list of names''', required=False, default='')
pars.add_argument('-ts', metavar='<str>', type=str, help='''list of names''', required=False, default='')
pars.add_argument('-reads_length', metavar='<int>', type=int, help='''reads length''', required=False, default=150)
pars.add_argument('-py', metavar='<str>', type=str, help='''name of python''', required=False, default='python3')


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Unsupported value encountered.')


pars.add_argument('-scaffold', type=str2bool, nargs='?', const=True, help='''make scaffold''', default=True)
pars.add_argument('-fast', type=str2bool, nargs='?', const=True, help='Turn on or turn off flag', default=True)
pars.add_argument('-force', type=str2bool, nargs='?', const=True, help='force do sth', default=False)
args = pars.parse_args()


# 判断格式，读取文件的函数
def judge_type(path):
    suffix = os.path.splitext(path)[-1].lower()
    suffix_dict = {'.gz': 0, '.fq': 1, '.fastq': 1, 'fa': 2, '.fas': 2, '.fasta': 2}
    return suffix_dict.get(os.path.splitext(path)[-1].lower(), 3)


def readrec(onefile):
    if judge_type(onefile) == 0:
        return SeqIO.parse(gzip.open(onefile, "rt"), "fastq")  # 读取二代测序数据的每一条reads
    elif judge_type(onefile) == 2:
        return SeqIO.parse(open(onefile), "fasta")  # fasta文件的每一条序列
    elif judge_type(onefile) == 1:
        return SeqIO.parse(open(onefile), "fastq")  # 读取二代测序数据的每一条reads
    else:
        print("Error: your input may invalid, please check your input.")
        sys.exit(1)


# 创建哈西字典,将reference存入字典
def gethashdict_v4(fileslist, merSize, kmer_dict, get_rc=False, clear=True, pos=False):
    gene_num = 0
    for file in fileslist:
        if os.path.isfile(file) == False: continue
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        file_name = os.path.basename(file)
        file_name = file_name.split('.')[0]
        if clear:
            with open(os.path.join(args.o, file_name + ".fasta"), "w") as f: pass
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            gene_num += 1
            refseq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            for j in range(0, len(refseq) - merSize + 1):
                temp_list, kmer = [], refseq[j:j + merSize]
                if kmer in kmer_dict:
                    temp_list = kmer_dict[kmer]
                    temp_list[0] += 1
                    if pos: temp_list[1] += (j + 1) / len(refseq)
                    if file_name not in temp_list: temp_list.append(file_name)
                else:
                    temp_list = [1, (j + 1) / len(refseq), file_name] if pos else [1, 0, file_name]
                    kmer_dict[sys.intern(kmer)] = temp_list
            if get_rc:
                try:
                    refseq = assemble.reverse_complement_limit(refseq)
                except:
                    print(file, refseq)
                for j in range(0, len(refseq) - merSize + 1):
                    temp_list, kmer = [], refseq[j:j + merSize]
                    if kmer in kmer_dict:
                        temp_list = kmer_dict[kmer]
                        temp_list[0] += 1
                        if pos: temp_list[1] -= (j + 1) / len(refseq)
                        if file_name not in temp_list: temp_list.append(file_name)
                    else:
                        temp_list = [1, - (j + 1) / len(refseq), file_name] if pos else [1, 0, file_name]
                        kmer_dict[sys.intern(kmer)] = temp_list
            if args.t > 0: print('Mem.:', round(sys.getsizeof(kmer_dict) / 300 / 500 / 1024, 2), 'G, Num. of Seq:',
                                 gene_num, end='\r')
        infile.close()


def get_ref_info(fileslist):
    global ref_length_dict
    for file in fileslist:
        ref_seq_count = 0
        if judge_type(file) == 2:
            file_name = os.path.basename(file)
            file_name = file_name.split('.')[0]
            ref_path_dict[file_name] = file
            for rec in readrec(file):
                refseq = ''.join(filter(str.isalnum, str(rec.seq).upper()))
                ref_seq_count += 1
                ref_length_dict[file_name] += len(refseq)
            ref_length_dict[file_name] = int(ref_length_dict[file_name] / max(ref_seq_count, 1))


def filter_reads(_dict, _kmer, stepSize, readseq, ignore_rc=False, single_file=False):
    tmplist = []
    for j in range(0, len(readseq) - _kmer + 1 - stepSize, stepSize):
        kmer = readseq[j:j + _kmer]
        if kmer in _dict:
            tmplist.extend(_dict.get(kmer)[2:])
            if single_file: break
    if not ignore_rc:
        temp_seq = assemble.reverse_complement_limit(readseq)
        for j in range(0, len(temp_seq) - _kmer + 1 - stepSize, stepSize):
            kmer = temp_seq[j:j + _kmer]
            if kmer in _dict:
                tmplist.extend(_dict.get(kmer)[2:])
                if single_file: break
    return set(tmplist)


def bytes_str(input, is_bytes_type):
    return input.decode('utf-8') if is_bytes_type else input


def do_reads_filter(_dict, _kmer, file, _re_dic, t_count, ignore_rc=False):
    t1, t2, reads_count, setpSize = time.time(), 0, 0, args.s
    single_file = len(_re_dic) == 1
    bytes_type = file[-3:].lower() == ".gz"
    infile = gzip.open(file, 'r') if bytes_type else open(file, 'r')
    for _ in infile:
        reads_count += 1
        temp_rec = [_, infile.readline(), infile.readline(), infile.readline()]
        for file_name in filter_reads(_dict, _kmer, setpSize, bytes_str(temp_rec[1], bytes_type), ignore_rc,
                                      single_file):
            with open(os.path.join(args.o, file_name + ".fasta"), "a+") as outfile:
                if file_name in _re_dic:
                    _re_dic[file_name] += 2
                else:
                    _re_dic[file_name] = 1
                outfile.writelines(['>', bytes_str(temp_rec[0], bytes_type), bytes_str(temp_rec[1], bytes_type)])
        if reads_count * t_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            print('handled\t', reads_count * t_count // 1000000, 'm reads, ', round(t2, 2), 's/m reads', sep="",
                  end='\r')
    infile.close()


def do_pair_reads_filter_v3(_dict, _kmer, file_1, file_2, _re_dic, t_id, t_count, ignore_rc=False):
    t1, t2, reads_count, setpSize = time.time(), 0, 0, args.s
    single_file = len(_re_dic) == 1
    bytes_type = file_1[-3:].lower() == ".gz"
    if file_1 == file_2:
        print("Error: the paired reads file is same!")
        exit(1)
    infile_1 = gzip.open(file_1, 'r') if bytes_type else open(file_1, 'r')
    infile_2 = gzip.open(file_2, 'r') if bytes_type else open(file_2, 'r')
    for _ in infile_1:
        reads_count += 1
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
        if reads_count % t_count == t_id:
            for file_name in filter_reads(_dict, _kmer, setpSize, bytes_str(temp_rec1[1], bytes_type), ignore_rc,
                                          single_file):
                with open(os.path.join(args.o, file_name + ".fasta"), "a+") as outfile:
                    if file_name in _re_dic:
                        _re_dic[file_name] += 2
                    else:
                        _re_dic[file_name] = 1
                    outfile.writelines(['>', bytes_str(temp_rec1[0], bytes_type), bytes_str(temp_rec1[1], bytes_type)])
                    outfile.writelines(['>', bytes_str(temp_rec2[0], bytes_type), bytes_str(temp_rec2[1], bytes_type)])
        if reads_count * 2 % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            print('handled\t', reads_count * 2 // 1000000, 'm reads, ', round(t2, 2), 's/m reads', " " * 4, sep="",
                  end='\r')
    infile_1.close()
    infile_2.close()


def filter_ref(fileslist, name_list):
    for file in fileslist:
        if judge_type(file) == 2:
            file_name = os.path.basename(file).split('.')[0]
            with open(os.path.join(args.o, file_name + ".fasta"), 'w') as outfile:
                for rec in readrec(file):
                    if rec.name.split("_")[1].upper() in name_list:
                        SeqIO.write(rec, outfile, "fasta")


def mylog(*myargs):
    _file = os.path.join(args.o, 'log.csv')
    with open(_file, 'a') as out:
        for i in myargs[:-1]: out.write(str(i) + ',')
        out.write(str(myargs[-1]) + '\n')


def get_FileSize(filePath):
    fsize = os.path.getsize(filePath)
    fsize = fsize / float(1024 * 1024)
    return round(fsize, 2)


def Execute_subproc(cmd):
    ex = subprocess.Popen(cmd)
    out, err = ex.communicate()
    status = ex.wait()
    return 1


if __name__ == '__main__':
    # 初始化操作
    if not os.path.isdir(args.o): os.mkdir(args.o)
    reads_length = args.reads_length
    ref_length_dict = defaultdict(int)
    ref_path_dict = defaultdict(str)
    used_dict = {}
    timestamp = args.ts if args.t < 0 else str(time.time())
    # py3.8 对于 macOS,默认启动方式是spawn,m1上会出错
    if sys.platform in ('darwin'):
        multiprocessing.set_start_method('fork')
    manager = multiprocessing.Manager()
    ref_readscount_dic = manager.dict()
    t0 = time.time()
    # 初始化参考序列信息
    if args.mode in [0, 1, 2, 3, 4, 5]:
        if args.t > 0: print("Getting information from references...")
        if os.path.isdir(args.f):
            path_list = os.listdir(args.f)
            for i in range(len(path_list)):
                path_list[i] = os.path.join(args.f, path_list[i])
            get_ref_info(path_list)
        elif os.path.isfile(args.f):
            get_ref_info([args.f])
        else:
            pass

    if args.mode == 1 or args.mode == 0:
        # 构建哈西字典
        path_list = []
        if args.t > 0: print('======================== Filter =========================')
        if os.path.isdir(args.f):
            for i in os.listdir(args.f):
                if judge_type(i) == 2:
                    path_list.append(os.path.join(args.f, i))
            gethashdict_v4(path_list, args.k1, used_dict, args.fast)
        elif os.path.isfile(args.f):
            gethashdict_v4([args.f], args.k1, used_dict, args.fast)
        else:
            print('Reference is invalid. Please check the path!')
            sys.exit(1)
        print("\nHash dictionary has been made")
        # print(len(used_dict))
        # 过滤-非配对文件
        if args.q:
            if args.t0 == 1:
                for t in range(len(args.q)):
                    do_reads_filter(used_dict, args.k1, args.q[t], ref_readscount_dic, len(args.q), args.fast)
            else:
                thread_list = []
                for t in range(len(args.q)):
                    p = Process(target=do_reads_filter,
                                args=(used_dict, args.k1, args.q[t], ref_readscount_dic, len(args.q), args.fast))
                    thread_list.append(p)
                for t in thread_list:
                    t.start()
                for t in thread_list:
                    t.join()
        # 过滤-配对
        # 硬盘速度瓶颈，限制2线程
        if args.q1 and args.q2:
            thread_list = []
            for i in range(len(args.q1)):
                if args.t0 == 1:
                    do_pair_reads_filter_v3(used_dict, args.k1, args.q1[i], args.q2[i], ref_readscount_dic, 0, 1,
                                            args.fast)
                else:
                    for t in (range(args.t0)):
                        p = Process(target=do_pair_reads_filter_v3, args=(
                        used_dict, args.k1, args.q1[i], args.q2[i], ref_readscount_dic, t, args.t0, args.fast))
                        thread_list.append(p)
                    for t in thread_list:
                        t.start()
                    for t in thread_list:
                        t.join()
        filter_count = 0
        for key, value in ref_readscount_dic.items():
            filter_count += 1
            # print('filter', key, filter_count,'/',len(ref_length_dict)," " * 16)
            mylog('filter', key, value, ref_length_dict[key], round(value / ref_length_dict[key] * reads_length, 2),
                  args.k1)
        t1 = time.time()
        print('Time used:', t1 - t0, " " * 32)
    # 超过512x重新梯度kmer过滤
    if args.mode == 2 or args.mode == 0:
        for i in range(1, args.t):
            command = [args.py, "miner.py", "-f", args.f, "-k1", args.k1, "-k2", args.k2, "-o", args.o, "-t", -i,
                       "-mode", "2",
                       "-limit_count", args.limit_count,
                       "-limit_min_length", args.limit_min_length,
                       "-limit_max_length", args.limit_max_length,
                       "-ts", timestamp,
                       "-reads_length", args.reads_length]
            _thread_subproc = threading.Thread(target=Execute_subproc, args=(list(map(str, command)),))
            _thread_subproc.start()

        args.s = 1
        if not os.path.isdir(os.path.join(args.o, 'big_reads')): os.mkdir(os.path.join(args.o, 'big_reads'))
        for key in ref_length_dict:
            if os.path.isfile(os.path.join(args.o, key + ".fasta")):
                if (args.mode == 2 or ref_readscount_dic.get(key, 0) / ref_length_dict[
                    key] * reads_length > 512) and get_FileSize(os.path.join(args.o, key + ".fasta")) > 8:
                    shutil.move(os.path.join(args.o, key + ".fasta"), os.path.join(args.o, 'big_reads', key + ".fasta"))
        filter_count = 0
        for key in ref_length_dict:
            filter_count += 1
            big_reads_path = os.path.join(args.o, 'big_reads', key + ".fasta")
            if os.path.isfile(big_reads_path):
                if os.path.isfile(os.path.join(args.o, key + ".fasta")) == False:
                    add_k = 2
                    while True:
                        print("re-filter:", ref_path_dict[key], 'with k =', args.k1 + add_k)
                        ref_readscount_dic[key] = 0
                        used_dict = {}
                        gethashdict_v4([ref_path_dict[key]], args.k1 + add_k, used_dict, True, True)
                        do_reads_filter(used_dict, args.k1 + add_k, big_reads_path, ref_readscount_dic, 1, True)
                        mylog('filter', key, ref_readscount_dic[key], ref_length_dict[key],
                              round(ref_readscount_dic[key] / ref_length_dict[key] * reads_length, 2), args.k1 + add_k)
                        if ref_readscount_dic[key] / ref_length_dict[key] * reads_length < 512 or get_FileSize(
                                os.path.join(args.o, key + ".fasta")) < 8 or args.k1 + add_k >= 127:
                            break
                        else:
                            add_k += 2
        # 等待其他进程完结
        if args.t > 0:
            wait_results = True
            while wait_results:
                wait_results = False
                for i in range(1, args.t):
                    if os.path.isfile(timestamp + str(-i)) == False:
                        wait_results = True
                        time.sleep(1)
                        break
            for i in range(1, args.t):
                os.remove(timestamp + str(-i))
        else:
            with open(timestamp + str(args.t), "w") as outfile:
                pass

    if args.mode == 3 or args.mode == 0:
        if args.t > 0:
            if not os.path.isdir(os.path.join(args.o, 'contig')):
                os.mkdir(os.path.join(args.o, 'contig'))
            if not os.path.isdir(os.path.join(args.o, 'short_contig')):
                os.mkdir(os.path.join(args.o, 'short_contig'))
            else:
                shutil.rmtree(os.path.join(args.o, 'short_contig'))
                os.mkdir(os.path.join(args.o, 'short_contig'))
            if not os.path.isdir(os.path.join(args.o, 'scaffold')):
                os.mkdir(os.path.join(args.o, 'scaffold'))
            else:
                shutil.rmtree(os.path.join(args.o, 'scaffold'))
                os.mkdir(os.path.join(args.o, 'scaffold'))

        for i in range(1, args.t):
            command = [args.py, "miner.py", "-f", args.f, "-k1", args.k1, "-k2", args.k2, "-o", args.o, "-t", -i,
                       "-mode", "3",
                       "-change_seed", args.change_seed,
                       "-limit_count", args.limit_count,
                       "-limit_min_length", args.limit_min_length,
                       "-limit_max_length", args.limit_max_length,
                       "-scaffold", args.scaffold,
                       "-ts", timestamp,
                       "-reads_length", args.reads_length]
            _thread_subproc = threading.Thread(target=Execute_subproc, args=(list(map(str, command)),))
            _thread_subproc.start()
        if args.t > 0: print('======================== Assemble =========================')
        used_dict = {}
        gc.collect()
        assemble_count, failed_count = 0, 0
        for key, value in ref_length_dict.items():
            limit = args.limit_count
            depth = 0
            assemble_count += 1
            contig_path = os.path.join(args.o, "contig", key + ".contig.fasta")
            short_contig_path = os.path.join(args.o, "short_contig", key + ".contig.fasta")
            if os.path.isfile(contig_path) == False and os.path.isfile(short_contig_path) == False:
                with open(contig_path, 'w') as out:
                    pass
                write_contig = False
                write_scaffold = args.scaffold
                print('assemble', key, assemble_count, '/', len(ref_length_dict), " " * 32)
                if os.path.isfile(os.path.join(args.o, key + ".fasta")) == False:
                    if os.path.isfile(contig_path): os.remove(contig_path)
                    continue
                # 获取种子列表
                ref_dict, filtered_dict = {}, {}
                gethashdict_v4([ref_path_dict[key]], args.k2, ref_dict, True, False, True)
                kmer_sum = assemble.make_assemble_hashdict([os.path.join(args.o, key + ".fasta")], args.k2,
                                                           filtered_dict, ref_dict)
                # 不在过滤出的reads里的置0
                for i in ref_dict:
                    if i not in filtered_dict:
                        ref_dict[i][0] = 0
                # limit为-1时根据count的分布计算新的limit
                if limit < 0:
                    if limit < 0:
                        limit, depth = assemble.dynamic_limit_v2(filtered_dict, value, kmer_sum, 2)
                        if depth > 20: limit = max(limit, assemble.dynamic_limit_v1(filtered_dict, 8, 64))
                # 保留与参考序列一致的kmer，不受limit限制
                if limit > 0:
                    _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit]
                    for x in _filter:
                        if filtered_dict[x][1][0] == 0:
                            del filtered_dict[x]
                    _filter = []
                if not filtered_dict:
                    failed_count += 1
                    if os.path.isfile(contig_path): os.remove(contig_path)
                    with open(short_contig_path, 'w') as out:
                        pass
                    print('Could not get enough reads from filter. Estimate depth:', depth, ' ' * 16)
                    continue
                # 根据序列的相似度计算权重
                cur_weight = assemble.dynamic_weight(filtered_dict)
                seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if
                             ref_dict[x][0] > 0 and ref_dict[x][1] > 0]
                list.sort(seed_list, key=lambda x: x[1], reverse=True)

                if not seed_list:
                    failed_count += 1
                    if os.path.isfile(contig_path): os.remove(contig_path)
                    with open(short_contig_path, 'w') as out:
                        pass
                    print('Could not get enough seeds. Estimate depth:', depth, " " * 16)
                    continue
                best_seed = seed_list[0]
                contig = assemble.get_contig_v5(filtered_dict, best_seed[0], weight=cur_weight)[0]
                best_contig = contig
                mylog('assemble', key, ref_readscount_dic.get(key, 0), value, len(contig), best_seed[0], args.k2, limit,
                      depth)
                # 如果序列太短
                if len(contig) / value < args.limit_min_length:
                    # 更换seed
                    org_file = os.path.join(args.o, key + ".fasta")
                    change_times, best_len = 0, len(contig)
                    last_seed = seed_list[0][0]
                    for seed in seed_list[1:]:
                        if last_seed[0][1:] == seed[0][:-1] or last_seed[0][2:] == seed[0][:-2]:
                            last_seed = seed
                            continue
                        last_seed = seed
                        print("Use another seed:", seed[0][:8], '...', change_times, " " * 16, end='\r')
                        contig = assemble.get_contig_v5(filtered_dict, seed[0], weight=cur_weight)[0]
                        if len(contig) > best_len: best_len, best_seed = len(contig), seed
                        mylog('assemble', key, ref_readscount_dic.get(key, 0), value, len(contig), seed[0], args.k2,
                              limit, depth)
                        change_times += 1
                        if len(best_contig) < len(contig): best_contig = contig
                        if len(best_contig) / value >= args.limit_min_length:
                            write_contig = True
                            break
                        if change_times >= int(args.change_seed): break
                    if len(best_contig) / value < args.limit_min_length:
                        failed_count += 1

                        if os.path.isfile(contig_path): os.remove(contig_path)
                    else:
                        write_contig = True
                # 如果序列太长
                elif len(contig) / value > args.limit_max_length:
                    # print(key, "Warning: Too long contig. Ref length:", value, ", contig length:", len(contig))
                    # 加大limit
                    for temp_limit in range(limit + 1, 33):
                        # 对filtered_dict再次精简
                        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= temp_limit]
                        for x in _filter:
                            if filtered_dict[x][1][0] == -1:
                                del filtered_dict[x]
                        _filter = []
                        # 此处更新权重影响很小，不更新
                        # cur_weight = assemble.dynamic_weight(filtered_dict)
                        contig = assemble.get_contig_v5(filtered_dict, best_seed[0], weight=cur_weight)[0]
                        mylog('assemble', key, ref_readscount_dic.get(key, 0), value, len(contig), best_seed[0],
                              args.k2, temp_limit, depth)
                        best_contig = contig
                        if len(contig) / value < args.limit_max_length:
                            break
                    write_contig = True
                else:
                    write_contig = True
                if write_scaffold:
                    seed_list.insert(0, best_seed)
                    scaffold = assemble.get_scaffold(filtered_dict, seed_list, limit, value, weight=cur_weight)
                    with open(os.path.join(args.o, 'scaffold', key + ".scaffold.fasta"), 'w') as out:
                        scaffold_length = len(scaffold.replace("N", ""))
                        out.write(
                            '>' + os.path.split(args.o)[-1] + '_scaffold_' + str(len(scaffold.replace("N", ""))) + '\n')
                        out.write(scaffold + '\n')
                        if not write_contig:
                            if scaffold_length / value > args.limit_min_length:
                                # print(key, " Try scaffold. Ref length:", value, ", best contig length: ", len(best_contig), ", scaffold length: ", scaffold_length, ", depth: ", depth)
                                failed_count -= 1
                            else:
                                print(key, " Failed. Ref length: ", value, ", best contig length: ", len(best_contig),
                                      ", scaffold length: ", scaffold_length, ", depth: ", depth, sep="")
                if write_contig:
                    with open(contig_path, 'w') as out:
                        out.write('>' + os.path.split(args.o)[-1] + '_contig_' + str(len(best_contig)) + '\n')
                        out.write(best_contig + '\n')
                else:
                    with open(short_contig_path, 'w') as out:
                        if not write_scaffold:
                            print(key, " Failed. Ref length: ", value, ", best contig length: ", len(best_contig),
                                  ", depth: ", depth, sep='')
                        out.write('>' + os.path.split(args.o)[-1] + '_short_contig_' + str(len(best_contig)) + '\n')
                        out.write(best_contig + '\n')
                    if os.path.isfile(contig_path): os.remove(contig_path)
                ref_dict, filtered_dict = {}, {}
                gc.collect()
        # 合并线程计数

        if args.t > 0:
            sum_failed_count = failed_count
            wait_results = True
            while wait_results:
                wait_results = False
                for i in range(1, args.t):
                    if os.path.isfile(timestamp + str(-i)) == False:
                        wait_results = True
                        time.sleep(1)
                        break
            for i in range(1, args.t):
                with open(timestamp + str(-i), "r") as infile:
                    sum_failed_count += int(infile.readline())
                os.remove(timestamp + str(-i))
            print("Assemble end. ", sum_failed_count, '/', assemble_count, 'failed', " " * 32)
        else:
            with open(timestamp + str(args.t), "w") as outfile:
                outfile.write(str(failed_count))

    if args.mode == 4:
        results_list = []
        muscle_exe = r"./muscle5"
        with open(args.list, "r") as infile:
            results_list = infile.readlines()
        Combine_count = 0
        for key, value in ref_length_dict.items():
            Combine_count += 1
            print('Combine', key, Combine_count, '/', len(ref_length_dict), " " * 16)
            in_file = os.path.join(os.path.join(args.o, key + ".combined.fasta"))
            out_file = os.path.join(os.path.join(args.o, key + ".aligned.fasta"))
            missed = False
            with open(os.path.join(args.o, key + ".combined.fasta"), 'w') as outfile:
                for result in results_list:
                    result = result.strip()
                    contig_file = os.path.join(result, "contig", key + ".contig.fasta")
                    scafflod_file = os.path.join(result, 'scaffold', key + ".scaffold.fasta")
                    if os.path.isfile(contig_file):
                        with open(contig_file, 'r') as infile:
                            outfile.writelines((infile.readlines()))
                    elif os.path.isfile(scafflod_file):
                        with open(scafflod_file, 'r') as infile:
                            # outfile.write('>'+result+'_scaffold\n')
                            # infile.readline()
                            # outfile.write(infile.readline())
                            outfile.writelines((infile.readlines()))
                    else:
                        missed = True
                        print('missing', contig_file)
                for rec in readrec(ref_path_dict[key]):
                    refseq = ''.join(filter(str.isalnum, str(rec.seq).upper()))
                    if refseq: outfile.write('>' + rec.name + '\n' + refseq + '\n')
            if missed == False or args.force == True:
                print(key, 'is aligning using muscle ...', " " * 16)
                muscle_cline = muscle_exe + " -align " + in_file + " -output " + out_file
                os.system(muscle_cline)
        t1 = time.time()
        print('Time used:', t1 - t0, " " * 32)

    if args.mode == 5:
        results_list = []
        trimal_exe = r"./trimal"
        Combine_count = 0
        for key, value in ref_length_dict.items():
            Combine_count += 1
            print('Trim', key, Combine_count, '/', len(ref_length_dict), " " * 16)
            in_file = os.path.join(os.path.join(args.o, key + ".aligned.fasta"))
            out_file = os.path.join(os.path.join(args.o, key + ".trimed.fasta"))
            if os.path.isfile(in_file):
                print(key, 'is triming using trimal ...', " " * 16, end='\r')
                trimal_cline = trimal_exe + " -in " + in_file + " -out " + out_file + " -gappyout"
                os.system(trimal_cline)
        t1 = time.time()
        print('Time used:', t1 - t0, " " * 32)

    if args.mode == 6:
        name_list = []
        print("Making new references...")
        with open(args.list, "r") as infile:
            for i in infile:
                name_list.append(i.strip().upper())
        if os.path.isdir(args.f):
            path_list = os.listdir(args.f)
            for i in range(len(path_list)):
                path_list[i] = os.path.join(args.f, path_list[i])
            filter_ref(path_list, name_list)
        elif os.path.isfile(args.f):
            filter_ref([args.f], name_list)
        else:
            pass
        t1 = time.time()
        print('Time used:', t1 - t0, " " * 32)
    t1 = time.time()
    if args.t > 0: print('Tatal Time used:', t1 - t0, " " * 32)
