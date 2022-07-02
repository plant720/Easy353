# Easy353 —— Short Manual

<center>zzhen</center><br><center>zzhen0302@163.com</center><br><center><font size="4">May 22,2022</font></center>

[TOC]

## 1. About Easy353

Easy353 is a tool for recovering Angiosperms353 gene set(AGS), which can filter and de novo assemble reads from sequencing data based on data from https://treeoflife.kew.org/, helping users capture AGS accurately and effectively.

**When using Easy353 please cite:**

XXXXX

## 2. Easy353 and Graphical User Interface

Easy353 is an easy-to-use pure Python program and is designed to be compatible with versions higher than 3.6. 

**Easy353 has a command line version and a GUI version**; users on Windows, macOS, and Linux computers could run it with cmd,and we also provide a user-friendly graphical interface for Windows and macOS users.

## 3. Downloading Easy353

Easy353 is open source under **MIT**. It is distributed via plant720 github repository: https://github.com/plant720/Easy353 where you can always download the most up to date version. Make sure to watch the github repository to remain up to date regarding code changes.

Version numbers follow the notation **x.y.z** where **x** changes with major code reorganizations, **y** changes when new features are added and **z** changes with bug fixes.

## 4. Installing Easy353

Easy353 comes in two different flavors: Easy353-cmd and Easy353-GUI. Thus, before installation you need to know on what kind of system you intend to execute Easy353. The Easy353-cmd can be run at popular systems(Linux, macOS and Windows); while the easy353-GUI is designed for Windows and macOS.

When working with small amounts of data on a Windows or Mac computer, we advice using the Easy353-GUI; larger amounts of data should be processed using the Easy353-cmd.

### Get Easy353-GUI

You can download Easy353-GUI from https://github.com/plant720/Easy353/release and Easy353-GUI can be run directly by double-clicking.

### Easy353-cmd

There are several generally 2 ways to install Easy353-cmd:

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

## 5. Input File Format

### 5.1 Reads Files

Reads files即通过二代测序获得的以fastq格式存储的测序数据。

```txt
@ST-E00600:58:HKYNGALXX:1:1101:1434:1000 1:N:0:GAGTTCGA
ATTGGGCACGACACGAAACGAAATTTTTGTTAACAGGAATTTAGAAGTAGATATTCAAGTTATGGTTGGTTATTAAGGTTATCTTATCAAAAGATAGAGCAAGATGGTTATGCTTTTTAAGGCTAAGGAAAATGTCAGCATCTACTTCAT
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ77JJJJJ7J7JJJ7JJJJJJJJ7JJJJJJ7JJJJJJJJJJJ-JJJJJJ-JJJJJ-J7-JJ
```

### 5.2 Reference files

参考序列是以fasta格式存储DNA序列的文本文件。每个参考文件包含一个以上的DNA序列——来自近源物种的同源基因。

每个序列的ID应该包括物种名和基因名，用`_`字符分开。示例如下：

```
>Litchi_chinensis_6946
CAAACCAAAATACACAATATAGGGGCAACACTTGTTGGGGTTGATAAATTTGGTAACAAGTATTATGAGAAACTTGGAAGACATAGGTGGGTTGAATATGCAGAGAAAG
>Tristiropsis_sp._6946
GGAAGACATAGGTGGGTTGAATATGCAGAGAAAGGTCGCTACAATGCCTCTCAGGTGCCTCCGGAATGGCATGGTTGGCTTCACTTTATAACTGATCACACAGGAGATG
```

## 6. Parameters

The most 



### 6.1 Basic parameters

The following options need to be specified, in the following order:

* `fq_file_1` and `fq_file_2`: The input files with paired-end reads, given in FASTQ format. 
* `unpaired_fq_file`: The input file with unpaired-end reads.

* `reference`: The file or 
* `output_dir`: Output directory

### 6.2 General parameters

The following option 

* `filter_kmer`:
* `assemble_kmer`:
* `filter_thread`
* `assemble_thread`

### 6.3 Advanced parameters

* `step_length`
* `ref_number`
* `change_seed`
* `kmer_limit`
* `minimum_length_ratio`
* `maximum_length_ration`
* `reference_number`
* `fast`



## 7. Output

## Running the pipeline

### 程序文件构成

>  程序主要是由 main.py util.py filter.py assembly.py 等文件构成

* `main.py`:软件运行的主程序。通过`python main.py -h`可以查看程序的参数

  ![image-20220522171614490](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220627154911.png)

* `filter.py`:用于从测序数据过滤与目标基因相关的reads

* `assembly.oy`:用于组装过滤出的reads

* `utils.py`:存放`filter.py`和`assembly.py`通用的函数

### 运行示例

以拟南芥*Arabidopsis thaliana*的转录组测序作为测序数据，十字花科被子植物353基因的fasta文件作为参考文件。其中十字花科被子植物353基因从是从https://treeoflife.kew.org下载并处理得到的。

```shell
python easy353.py -1 SRR18391637.sra_1.fastq -2 SRR18391637.sra_2.fastq -r ref_Brassicaceae -o result -k1 31 -k2 41
```

```shell
# 参数解释
-1 -2: 二代测序的双端测序数据
-r: 参考序列文件。当需要捕获多个基因时，需要将不同基因的参考序列文件存放在文件夹中，并给出文件夹路径；当只需要筛选单个基因时，只需给出单个文件路径即可
-o: 指定输出文件夹名
-k1: 用于过滤测序数据的k值
-k2: 用于组装读长时的k值
```

### 输出文件

> 输出文件夹主要包括reads、contig、scaffold、big_reads 、short_contig等文件夹及assembly_log.txt、filter_log.txt及result_log.txt等文件

* contig:用于存放contig的文件夹。该文件夹中的不同文件是根据不同基因的参考文件过滤出的reads组装出的最优单一contig。当组装出的contig长度较长，可作为捕获出的目标基因使用。contig文件中的contig长度在目标基因平均长度的100%-200%。
* short_contig:当组装出的contig长度没有达到目标基因长度的50%，会被放入short_contig文件夹中。
* scaffold:存放scaffold文件，scaffold是使用filtered reads组装出的多条contig，根据在参考序列中的相对位置和方法连接得到的。
* reads:存放测序数据中与目标基因参考序列局部相似的reads。
* big_reads:当单一基因通过过滤得到的reads数量过多时，存入big_reads文件夹，通过扩大k值再次进行过滤。

> 一般来说，当拼接得到的contig数量较多、长度较长时，可以采用contig文件夹中的contig作为目标基因使用。当组装得到的contig长度较短，可以考虑采用scafflod中得到的结果进行使用，但是因为间隙（gap）的存在，需要进行修剪才可使用。

* filter_log.txt:记录从测序数据中过滤reads的信息，主要包括过滤基因名、过滤出的reads数量、目标基因平均长度、覆盖度和K值
* assembly_log.txt:存放组装reads过程中的信息，主要包括组装基因名、组装出的contig长度、目标基因平均长度、组装过程中的seed、K值
* result_lof.txt:存放程序运行的结果信息，主要包括基因名、目标基因长度、过滤出的reads数量、contig长度、scaffold长度等