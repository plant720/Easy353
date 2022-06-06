import argparse
import sys
import os
import subprocess
import glob
from Bio import AlignIO
import csv
import shutil


pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                               description=''' combine and align seqs ''')
group = pars.add_mutually_exclusive_group(required=True)
group.add_argument('-i', metavar='<str>', type=str, help='''input pathes.''', nargs="+")
group.add_argument('-l', metavar='<str>', type=str, help='''input pathes list:spname    sp_contig_path''')
pars.add_argument('-g', metavar='<str>', type=str, help='''input genelist.''', required=True)
pars.add_argument('-o', metavar='<str>', type=str, help='''out path.''', required=False, default="contig_deal")
pars.add_argument('-t', metavar='<int>', type=int, help='''t.''', required=False, default=24)
pars.add_argument('-s', help='''stats mafft result.''', action='store_true', default=False)
args = pars.parse_args()


def get_files(ref):
    file_path_list = []
    for root, dirs, files in os.walk(ref):
        for file in files:
            file_path = os.path.join(root, file)
            file_path_list.append(file_path)
    return file_path_list


def is_exist(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0
    elif os.path.isdir(file):
        files = get_files(file)
        if files == []:
            flag = 0
        else:
            flag = 1
            # for i in files:
            #     if os.path.getsize(i) > 0:
            #         continue
            #     else:
            #         flag = 0
            #         break
    else:
        flag = 0
    return flag


def get_basename(file):
    if is_exist(file):
        basename = os.path.basename(file)
        if ".fasta" in basename:
            basename = basename.split(".")[0]
        elif ".fas" in basename:
            basename = basename.split(".")[0]
        elif ".fa" in basename:
            basename = basename.split(".")[0]
        else:
            basename = basename
        return basename


def getabpath(pathlist):
    abpathlist = [os.path.abspath(i) for i in pathlist]
    return abpathlist


def find_file(path, pattern):
    if is_exist(path):
        file = glob.glob(os.path.join(path, f"*{pattern}*"))
        if file != [] and is_exist(file[0]):
            # if is_exist(file[0])
            return os.path.abspath(file[0])
        else:
            return False
    else:
        print("please chack your input path!")
        sys.exit(0)


def mkdir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
        outpath = os.path.abspath(path)
    else:
        os.makedirs(path, exist_ok=True)
        outpath = os.path.abspath(path)
    return outpath


def info_twofasta(seq1, seq2):
    allmatch = 0
    match = 0
    tmps = ""
    for i, j in zip(seq1, seq2):
        if i == "-" or j == "-":
            pass
        else:
            allmatch += 1
            if i == j:
                tmps += i
                match += 1
    if allmatch != 0 and match != 0:
        if tmps in seq1:
            insert = 0
        else:
            insert = 1
        identity = float(100 * match / allmatch)
        coverage1 = float(100 * allmatch / len(str(seq1.seq).replace('-', '')))
        coverage2 = float(100 * allmatch / len(str(seq2.seq).replace('-', '')))
        return identity, coverage1, coverage2, insert
    else:
        return False


def mylog(path, sth):
    with open(path, "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(sth)


def getcmdreturn(cmd):
    result = subprocess.getoutput(cmd)
    return result


def parse_align(alignfiles, outpath):
    for file in alignfiles:
        align = AlignIO.read(file, "fasta")
        name_list = []
        count = 1
        for i in align:
            for j in align[count:]:
                temp = {i.name, j.name, }
                if temp not in name_list:
                    name_list.append(temp)
                else:
                    continue
                a = info_twofasta(i, j)
                if a:
                    info = [get_basename(file), i.name, j.name, a[0], a[1], a[2], a[3]]
                    path = os.path.join(outpath, "mafft.xls")
                    mylog(path, info)
                    # print(i.name,j.name,a[0],a[1],a[2],sep="\t")
            count += 1


def combine(gene_list, path, inputpath=args.i, pathlist=args.l):
    outpath = mkdir(path)
    genefile = os.path.abspath(args.g)
    combinecout = 0
    if inputpath:
        inputabpath = getabpath(inputpath)

        for gene in gene_list:

            combine_file = os.path.join(outpath, gene + ".combined.fasta")
            mafft_file = os.path.join(outpath, gene + ".aligned.fasta")
            runfile = os.path.join(outpath, "run.sh")
            # print(runfile)
            missed = False
            with open(combine_file, 'w') as comf:
                for inpath in inputabpath:
                    # print(inpath,gene)
                    contig_file = find_file(inpath, gene)
                    if contig_file != False:
                        with open(contig_file, 'r') as infile:
                            inlines = infile.readlines()
                            # inlines[0] = "_".join(inlines[0].split("_")[:2]) + "\n"
                            comf.writelines(inlines)
                    else:
                        missed = True
                        os.remove(combine_file)
                        break
            # print(missed)
            if missed == False:
                combinecout += 1
                with open(runfile, 'a') as runf:
                    command = f"mafft --maxiterate 1000 --localpair {combine_file} >{mafft_file}\n"
                    runf.writelines(command)
        if combinecout == 0:
            print("err,please chack your input path!")
            sys.exit(0)
        else:
            os.chdir(outpath)
            if subprocess.call(f"nohup ParaFly -c run.sh -CPU {args.t} ", shell=True) == 0:
                print("combine and mafft has done!")
    elif pathlist:
        pathdict = {}
        with open(pathlist, "r") as f:
            for line in f:
                L = line.strip().split(sep=None)
                # print(L)
                pathdict[L[0]] = os.path.abspath(L[1])
                # print(pathdict)
        for gene in gene_list:
            combine_file = os.path.join(outpath, gene + ".combined.fasta")
            mafft_file = os.path.join(outpath, gene + ".aligned.fasta")
            runfile = os.path.join(outpath, "run.sh")
            missed = False
            with open(combine_file, 'w') as comf:
                for sp, sppath in pathdict.items():
                    # print(sp,gene)
                    contig_file = find_file(sppath, gene)
                    if contig_file != False:
                        with open(contig_file, 'r') as infile:
                            inlines = infile.readlines()

                            inlines[0] = f">{sp}_{gene}" + "\n"
                            comf.writelines(inlines)
                    else:
                        missed = True
                        os.remove(combine_file)
                        break
            if missed == False:
                with open(runfile, 'a') as runf:
                    command = f"mafft --maxiterate 1000 --localpair {combine_file} >{mafft_file}\n"
                    runf.writelines(command)
        os.chdir(outpath)
        if subprocess.call(f"nohup ParaFly -c run.sh -CPU {args.t} ", shell=True) == 0:
            print("combine and mafft has done!")


if __name__ == '__main__':
    with open(args.g, "r") as infile:
        gene_list = [i.strip() for i in infile.readlines()]

    # 01.combine_mafft
    alloutpath = mkdir(args.o)
    com_mafft_path = os.path.join(alloutpath, "01.combine_mafft")
    combine(gene_list, com_mafft_path)
    mafftfiles = glob.glob(os.path.join(com_mafft_path, f"*aligned.fasta"))
    if args.s:
        parse_align(mafftfiles, com_mafft_path)
        mafftlog = os.path.join(com_mafft_path, 'mafft.xls')
        mafftgood = getcmdreturn(f"awk '$4==100 && $6==100 && $7==0' {mafftlog} |wc -l")
        print(mafftgood)
        # with open("/data/user/gyl/AT_dataset/get_353/test.result.xls","a") as f:
        # f.write(f"{args.o}\t{mafftgood}\n")
    os.chdir(alloutpath)
    # python compare.py -i result/contig/ ref_gold/ -g 353gene.list -s -o stats
