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

Using the setup.py, **you should have Python library setuptools installed** (`sudo apt install -y python-setuptools` or `sudo yum install -y python-setuptools`or `pip install setuptools`).

## 5. Input File Format

### 5.1 Reads Files

Easy353 assembles any High-throughput Sequencing(HTS) reads, given in the **FASTQ** format. Paired or unpaired reads are OK. And Easy353 can direclty read files compressed with gzip. Compressed files should end with `.gz`.

```txt
@ST-E00600:58:HKYNGALXX:1:1101:1434:1000 1:N:0:GAGTTCGA
ATTGGGCACGACACGAAACGAAATTTTTGTTAACAGGAATTTAGAAGTAGATATTCAAGTTATGGTTGGTTATTAAGGTTATCTTATCAAAAGATAGAGCAAGATGGTTATGCTTTTTAAGGCTAAGGAAAATGTCAGCATCTACTTCAT
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ77JJJJJ7J7JJJ7JJJJJJJJ7JJJJJJ7JJJJJJJJJJJ-JJJJJJ-JJJJJ-J7-JJ
```

### 5.2 Reference files

DNA sequences are stored in reference files in **FASTA** format. Each reference file contains one or more DNA sequence,which is the homologous gene from closely related species. Here we provide an example: https://github.com/plant720/Easy353Test/tree/master/data

Each sequence's ID should contain the gene name and species name, separated by the `_`.

```
>Litchi_chinensis_6946
CAAACCAAAATACACAATATAGGGGCAACACTTGTTGGGGTTGATAAATTTGGTAACAAGTATTATGAGAAACTTGGAAGACATAGGTGGGTTGAATATGCAGAGAAAG
>Tristiropsis_sp._6946
GGAAGACATAGGTGGGTTGAATATGCAGAGAAAGGTCGCTACAATGCCTCTCAGGTGCCTCCGGAATGGCATGGTTGGCTTCACTTTATAACTGATCACACAGGAGATG
```

## 6. Parameters

### 6.1 Parameters to be specified

The following options need to be specified:

* `fq_file_1` and `fq_file_2`: The input files with paired-end reads, given in FASTQ format. 
* `unpaired_fq_file`: The input file with unpaired-end reads.

* `reference`:  A directory recording 353 gene sequences from closely related species, given in FASTA format. 
* `output_dir`: The output directory. 

### 6.2 Parameters with default value

The following options have default values and do not need to be specified in general.

* `filter_kmer`: K-mer length setting for filtering, default value is 31. The `filter_kmer` is the length of the shared DNA segment between the reads and the reference sequences.  Easy353 looks for a shared K-length DNA segment (K-mer) between a read and a reference to assess whether they are related. Consequently, when the value of  `filter kmer` is set to 31, it means to search for reads that have 31bp in common with  the reference sequence.
* `assemble_kmer`: K-mer length setting for assembly, default value is 41.  The `assemble_kmer` is the length of the nodes in the de Bruijn graph (DBG). In reads assembly, the filtered reads are divided into K-mers; the k-mers are employed as the nodes of the DBG. It strongly depends on the input dataset. For Illumina reads(150bp) with sufficient coverage (> 40x), we have good results with k = 41.
* `kmer_limit`: Limit of K-mer count, default value is 2. The`kmer_limit` is used to remove erroneous, low-abundance K-mers. This parameter also strongly depends on the dataset. It corresponds to the smallest amount of times a correct k-mer appears in the reads. A typical value is 2, which means only K-mers with at least 3 occurrences will be used for assembly.  If the dataset has high coverage, try larger values.
* `filter_thread`: Threads setting for reads filtering, defalut value is 1. This value is advised to ideally not exceed 4. And for Windows and macOS, the `filter thread` should be set to 1.
* `assemble_thread`: Threads setting for reads assembly, defalut value is 4. 
* `reference_number`: The number of the reference sequences used to build hash table, and the full reference sequence will be used by default. When one target gene has more than 100 homologous reference sequences, we advise setting this `reference_number` to 100 since the amount of sequences in the fasta file will effect the memory size needed by Easy353.
* `fast`: Whether to use fast mode. 

* `step_length`
* `change_seed`
* `minimum_length_ratio`: The minimum ratio of the assembled gene length to reference average length, default value is 1.0. 
* `maximum_length_ration`

## 7. Output

containing the `filtered_reads`  and the `target_genes` . The directory, `filtered_reads`, contains filtered reads associated with target 353 genes. And the directory, `target_genes`, is the most important result directory of Easy353, which stores 353 gene sequences recovered by easy353.