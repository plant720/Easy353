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
# get a local copy of easy353 source code
git clone https://github.com/plant720/Easy353.git
# install the code 
cd Easy353
python setup.py install --user
```

Using the setup.py, you should have Python library setuptools installed (`sudo apt install -y python-setuptools` or `sudo yum install -y python-setuptools`or `pip install setuptools`).

## Test

* Download the simulation data of [*Glycine max*](https://github.com/plant720/Easy353Test/tree/master/data):

```shell
wget https://github.com/plant720/Easy353Test/blob/master/data/Gmax_sim_100_fir.fastq.gz
wget https://github.com/plant720/Easy353Test/blob/master/data/Gmax_sim_100_sec.fastq.gz
```

### Easy353-cmd

* After installation of Easy353 and download the simulation data of *Glycine max*, please download the AGS of related species as the reference

```shell
# download the AGS according to taxonomy: Glycine_max is one species of Fabaceae family
# The ref can be found at ref/353gene
build_database.py -o 353_ref_Fabaceae -c Fabaceae -t 10 -exclude Glycine_max -generate 
```

* then do the recover of  Angiosperms353 gene set(AGS)

```shell
easy353.py -1 Gmax_sim_100_fir.fastq.gz -2 Gmax_sim_100_sec.fastq.gz -r 353_ref_Fabaceae -o test_package -k1 31 -k2 41 -t1 4 -t2 4 -fast
```

You are going to get the same result as [here](https://github.com/plant720/Easy353Test/tree/master/result/Gmax_100_result).

### Easy353-GUI

* For Easy353-GUI, you should set the `Paired fq file 1`, `Paired fq file 2`,`Taxonomy` and `Output dir` as shown in the image below

![image-20220629211901122](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220629211901.png)

## User manual

A more complete manual is here: https://github.com/plant720/Easy353/docs/manual.pdf

## Contact

Questions, suggestions, comments, etc? You can email the developer for specific support: zzhen0302@163.com
