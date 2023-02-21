# Easy353
If you use the Easy353 software, please cite our manuscript: Zhen Zhang, Pulin Xie, Yongling Guo, Wenbin Zhou, Enyan Liu, Yan Yu. **Easy353: A tool to get Angiosperms353 genes for phylogenomic research**. *Molecular Biology and Evolution*. msac261 (2022). https://doi.org/10.1093/molbev/msac261.

Please also cite the dependencies if used:

Baker W.J., Bailey P., Barber V., Barker A., Bellot S., Bishop D., Botigue L.R., Brewer G., Carruthers T., Clarkson J.J., Cook J., Cowan R.S., Dodsworth S., Epitawalage N., Francoso E., Gallego B., Johnson M., Kim J.T., Leempoel K., Maurin O., McGinnie C., Pokorny L., Roy S., Stone M., Toledo E., Wickett N.J., Zuntini A.R., Eiserhardt W.L., Kersey P.J., Leitch I.J. & Forest F. **A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life**. *Systematic Biology*. 71: 301–319. https://doi.org/10.1093/sysbio/syab035.

## Introduction

Easy353 is a tool for recovering Angiosperms353 gene set(AGS), which can filter reads from high throughput sequencing data such as RNASeq and genome skimming and capture AGS accurately and effectively with our optimized reference-guided assembler. 

## Installation

Easy353 is an easy-to-use program based on Python 3 (3.6 and above), which is distributed with two user interfaces: a full graphical user interface (Easy353-GUI) and a command-line interface (Easy353). Users who have laptops or desktop computers could use the GUI version for Easy353. 
We also provide a command line for computers with large memory and multi-core CPUs.

### Easy353 with GUI

For Windows and macOS users, it is recommended to get the ALL-IN-ONE graphical user interface version of Easy353 (Easy353-GUI) from https://github.com/plant720/Easy353/releases.

![image-20220629211901122](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220629211901.png)

### Command line version of Easy353

There are several generally 2 ways to install Easy353:

* Option 1 **Using the setup.py**
* Option 2 **In situ configuration**

#### Option 1. Using the setup.py

The most up-to-date version of Easy353 is available at our GitHub site (https://github.com/plant720/Easy353). Users should use git to download the entire Easy353 repository and install the Easy353 using the setup.py. On macOS, users could open a command terminal by Applications->Utilities->Terminal.

```shell
# get a local copy of the easy353 source code
git clone https://github.com/plant720/Easy353.git
# install the code 
cd Easy353
python3 setup.py install --user
```

Using the setup.py, you should have Python library setuptools installed (`sudo apt install -y python-setuptools` or `sudo yum install -y python-setuptools` or `pip3 install setuptools`).

For some Linux nad macOS systems, after above commands you still cannot execute `build_database.py` and `easy353.py` in a new terminal directly, meaning `~/.local/bin` was not added to the \$PATH, you have to manually add `~/.local/bin` :

```shell
# add ~/.local/bin to PATH
echo "PATH=~/.local/bin:\$PATH" >> ~/.bashrc
source /.bashrc
```

#### Option 2. In situ configuration

Use git to download the entire Easy353 repository.

```shell
# Supposing you are going to install it at ~/Applications
mkdir ~/Applications # create directories if not exist
cd ~/Applications
git clone https://github.com/plant720/Easy353.git
```

Use the following commands to make Easy353 scripts executable

```shell
chmod 755 Easy353/build_database.py
chmod 755 Easy353/easy353.py
```

Add Easy353 to the $PATH.

```shell
echo "export PATH=~/Applications/Easy353:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

At last, install python libraries biopython, psutil, requests, and beautifulsoup4 using pip or conda.

```shell
# install required libs
pip install biopython psutil requests beautifulsoup4
```

## Test
The test files can be found at https://github.com/plant720/Easy353Test

* Download the simulation data of [*Glycine max*](https://github.com/plant720/Easy353Test/tree/master/Glycine_max):

```shell
wget https://raw.githubusercontent.com/plant720/Easy353Test/master/Glycine_max/Gmax_sim_1.fastq.gz
wget https://raw.githubusercontent.com/plant720/Easy353Test/master/Glycine_max/Gmax_sim_2.fastq.gz
```

### Easy353-cmd

* After installation of Easy353 and downloading the simulation data of *Glycine max*, please download the AGS of related species as the reference

```shell
# Download the AGS data according to taxonomy: Glycine max is a species from Glycine genus in Fabaceae, so we download species from Fabaceae as the reference.
# The reference sequences are downloaded from Kew Tree of Life Explorer (https://treeoflife.kew.org), so keep your devices connected to the network. And if you the build_database.py, please cite: https://doi.org/10.1093/sysbio/syab035.
build_database.py -o 353_ref_Fabaceae -c Fabaceae -t 10 -exclude Glycine_max -generate 
# The final reference sequences can be found at 353_ref_Fabaceae/353gene after downloading

## Explanation of parameters
-o: the output directory
-c: the taxonomy of species that used as reference
-t: the thread used to download files
-exclude: exclude species that are not used as reference
-generate: generate a CSV that records the info of downloaded species
```

* Then do the recovery of Angiosperms353 gene set(AGS)

```shell
# use easy353.py to filter and assemble reads to get target genes
easy353.py -1 Gmax_sim_1.fastq.gz -2 Gmax_sim_2.fastq.gz -r 353_ref_Fabaceae/353gene -o test_package -k1 31 -k2 41 -t1 1 -t2 4 -reference_number 100

## Explanation of parameters
-1 -2: the input files with paired-end reads, given in FASTQ format. 
-r: the reference directory
-o: the output directory
-k1: K-mer length setting for filtering
-k2: K-mer length setting for assembly
-t1: the threads setting for reads filtering
-t2: the threads setting for reads assembly
-reference_number: the number of the reference sequences used
```

* Now, you can view the result of Easy353 in the output directory

The output directory contains the `filtered_reads` and the `target_genes`. 

* The directory, `filtered_reads`, contains filtered reads associated with target 353 genes. 
* And the directory, `target_genes`, is the most important result directory of Easy353, which stores 353 gene sequences recovered by Easy353. 
  * The files under `target_genes` are recovered successfully, and the `unrecovered_genes`  contains 353 genes that were not recovered successfully because the recovered genes are short in length.

### Easy353-GUI

Please use the Easy353-gui.app or Easy353-gui.exe from https://github.com/plant720/Easy353/releases, the `Easy353-gui.py` is only used for development.

* For Easy353-GUI, you should set the `Paired fq file 1`, `Paired fq file 2`, `Taxonomy` ,and `Output dir` as shown in the image below

![image-20220629211901122](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220629211901.png)

## User manual

A more complete manual is here: https://github.com/plant720/Easy353/blob/master/docs/manual.pdf

## Contact

Questions, suggestions, comments, etc? You can email the developer for specific support: zzhen0302@163.com
