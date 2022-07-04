# Easy353

## Introduction

Easy353 is a tool for recovering Angiosperms353 gene set(AGS), which can filter and de novo assemble reads from sequencing data based on data from https://treeoflife.kew.org/, helping users capture AGS accurately and effectively.

## Installation

Easy353 is an easy-to-use pure Python program and is designed to be compatible with versions higher than 3.6. Users on Windows, macOS, and Linux computers could run Easy353 directly from the command line. We also provide a user-friendly graphical interface for Windows and macOS. 

### Easy353 with GUI

It is advised to get Easy353-GUI from https://github.com/plant720/Easy353/release if you use Windows or macOS. Tha app file for macOS and the exe file for Windows, which be run by double-clicking.

### Easy353 for command line

There are several generally 2 ways to install Easy353:

* Option 1 **Using conda**

- Option 2 **Using the setup.py**

#### Option 1. Using conda



#### Option 2. Using the setup.py

You should use git to download the entire Easy353 repository and install the Easy353 using the setup.py.

```shell
# get a local copy of the easy353 source code
git clone https://github.com/plant720/Easy353.git
# install the code 
cd Easy353
python setup.py install --user
```

Using the setup.py, you should have Python library setuptools installed (`sudo apt install -y python-setuptools` or `sudo yum install -y python-setuptools`or `pip install setuptools`).

## Test

* Download the simulation data of [*Glycine max*](https://github.com/plant720/Easy353Test/tree/master/data):

```shell
wget https://github.com/plant720/Easy353Test/raw/master/data/Gmax_sim_1.fastq.gz
wget https://github.com/plant720/Easy353Test/raw/master/data/Gmax_sim_2.fastq.gz
```

### Easy353-cmd

* After installation of Easy353 and downloading the simulation data of *Glycine max*, please download the AGS of related species as the reference

```shell
# Download the AGS data according to taxonomy: Glycine max is a species from Glycine genus in Fabaceae, so we download species from Fabaceae as the reference.
# The reference sequences are downloaded from https://treeoflife.kew.org/,so keep your devices connected to the network. And if you the build_database.py, please cite: https://doi.org/10.1093/sysbio/syab035.
build_database.py -o 353_ref_Fabaceae -c Fabaceae -t 10 -exclude Glycine_max -generate 
# The final reference sequences can be found at 353_ref_Fabaceae/353gene after downloading

## Explanation of parameters
-o: the output directroy
-c: the taxonomy of species that used as reference
-t: the thread used to download files
-exclude: exclud species that are not used as reference
-generate: generate a csv that records the info of downloaded species
```

* Then do the recovery of Angiosperms353 gene set(AGS)

```shell
# use easy353.py to filter and assemble reads to get target genes
easy353.py -1 Gmax_sim_1.fastq.gz -2 Gmax_sim_2.fastq.gz -r 353_ref_Fabaceae/353gene -o test_package -k1 31 -k2 41 -t1 1 -t2 4 -reference_number 100

## Explanation of parameters
-1 -2: the input files with paired-end reads, given in FASTQ format. 
-r: the reference directory
-o: the output directiry
-k1: K-mer length setting for filtering
-k2: K-mer length setting for assembly
-t1: the threads setting for reads filtering
-t2: the threads setting for reads assembly
-reference_number: the number of the reference sequences used
```

* Now, you can view the result of Easy353 in output directory

The output directory contains the `filtered_reads` and the `target_genes` . 

* The directory, `filtered_reads`, contains filtered reads associated with target 353 genes. 
* And the directory, `target_genes`, is the most important result directory of Easy353, which stores 353 gene sequences recovered by Easy353. 
  * The files under `target_genes` are recovered successfully, and the `unrecovered_genes`  contains 353 genes that were not recovered successfully, because the recovered genes are shorter in length.

### Easy353-GUI

Please use the Easy353-gui.app or Easy353-gui.exe from https://github.com/plant720/Easy353/release, the `Easy353-gui.py` is only used for developing.

* For Easy353-GUI, you should set the `Paired fq file 1`, `Paired fq file 2`, `Taxonomy` and, `Output dir` as shown in the image below

![image-20220629211901122](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220629211901.png)

## User manual

A more complete manual is here: https://github.com/plant720/Easy353/docs/manual.pdf

## Contact

Questions, suggestions, comments, etc? You can email the developer for specific support: zzhen0302@163.com
