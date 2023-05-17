<p align="center">
   <img src="https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/202304272252769.svg" alt="logo22" width="200";" />
</p>


### Easy353: An Efficient Tool Designed to Recover the Angiosperms353 Gene Set

**Notice: The Easy353 has been updated to v2.0.1, which is faster and more accurate.**

Please cite the following manuscript if you use the Easy353:

Zhen Zhang, Pulin Xie, Yongling Guo, Wenbin Zhou, Enyan Liu, Yan Yu. **Easy353: A tool to get Angiosperms353 genes for phylogenomic research**. *Molecular Biology and Evolution*. msac261 (2022). https://doi.org/10.1093/molbev/msac261.

Additionally, please cite the dependencies if used:

Baker W.J., Bailey P., Barber V., Barker A., Bellot S., Bishop D., Botigue L.R., Brewer G., Carruthers T., Clarkson J.J., Cook J., Cowan R.S., Dodsworth S., Epitawalage N., Francoso E., Gallego B., Johnson M., Kim J.T., Leempoel K., Maurin O., McGinnie C., Pokorny L., Roy S., Stone M., Toledo E., Wickett N.J., Zuntini A.R., Eiserhardt W.L., Kersey P.J., Leitch I.J. & Forest F. **A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life**. *Systematic Biology*. 71: 301–319. https://doi.org/10.1093/sysbio/syab035.

Easy353 is a tool specifically designed to recover the Angiosperms353 gene set (AGS). It effectively filters AGS-related reads from high-throughput sequencing data, and accurately recovers AGS using its optimized reference-guided assembler. 

![image-20230427223819794](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/202304272238784.png)

To run Easy353, two things are required as input: the sequencing data (reads in FASTQ format) and a set of reference sequences (orthologous genes from other species). The target genes of species X are included in the output file. **Notice: Easy353 can recover not only AGS but also user-specified target genes (e.g., chloroplast genes, ITS). However, if the target gene is a gene like ITS, the user must provide their own reference sequences. Using the script `build_database.py`, users can download the AGS reference sequences while recovering AGS.**

## Download and Installation

Easy353 is a user-friendly tool developed in Python 3 (3.6 and above), offering two interfaces: a full graphical user interface (Easy353-GUI) and a command-line interface (Easy353).

### Easy353-GUI

For Windows and macOS users, it is recommended to download the ALL-IN-ONE graphical user interface version of Easy353 (Easy353-GUI) from [Easy353/releases](https://github.com/plant720/Easy353/releases).

### Easy353-cmd

There are several generally 2 ways to install Easy353-cmd:

* Option 1 **Using the setup.py**
* Option 2 **In situ configuration**

#### Option 1. Using the setup.py

The most up-to-date version of Easy353 is available at our GitHub site [plant720/Easy353](https://github.com/plant720/Easy353). Users should use git to download the entire Easy353 repository and install the Easy353 using the `setup.py`. 

1. Clone the Easy353 repository:

```shell
# get a local copy of the easy353 source code
git clone https://github.com/plant720/Easy353.git
```

2. Install the code:

```shell
 # install the code 
 cd Easy353
 python3 setup.py install --user
```

3. For some Linux and macOS systems, after executing the above commands, you may not be able to run `build_database.py` and `easy353.py` directly in a new terminal, indicating that `~/.local/bin` has not been added to the `$PATH`. In this case, you have to manually add `~/.local/bin`:

```shell
# add ~/.local/bin to PATH
echo "PATH=~/.local/bin:\$PATH" >> ~/.bashrc
source /.bashrc
```

#### Option 2. In situ configuration

Alternatively, you can use the following commands to download and configure Easy353 in situ:

1. Clone the Easy353 repository and create a directory for installation:

```shell
# Assuming you want to install it at ~/Applications
mkdir ~/Applications # create directories if not exist
cd ~/Applications
git clone https://github.com/plant720/Easy353.git
```

2. Make the Easy353 scripts executable:

```shell
chmod +x Easy353/build_database.py
chmod +X Easy353/easy353.py
```

3. Add Easy353 to the `$PATH`.

```shell
echo "export PATH=~/Applications/Easy353:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

4. Install the required Python libraries, including biopython, psutil, requests, and beautifulsoup4:

```shell
# install required libs
pip install biopython psutil requests beautifulsoup4
```

## Usage

**Note: The detailed usage and tutorials can be found at [Easy353/wiki](https://github.com/plant720/Easy353/wiki)**.

After installation, Easy353 can be used to filter and assemble reads to recover target genes.

1. Prepare the input:

**To run Easy353, two things are required as the input: the sequencing reads in FASTQ format and a set of reference sequences, i.e. AGS.** The quality of the reference sequences strongly influences the Easy353's result. Hence, it is crucial to generate the reference sequences carefully. If the target sequence is AGS, users can use the `build_database.py` script to generate. However, if the target gene is other marker genes, the users need to generate the reference sequences independently. Here is an example of how to use our script to generate the reference:

```shell
# To ensure that the downloaded AGS data consists of sequences from closely related species, please download the sequences according to taxonomy.
# The reference sequences are downloaded from Kew Tree of Life Explorer (https://treeoflife.kew.org), so ensure your device is connected to the network.
build_database.py -o 353_ref_Fabaceae -c Fabaceae -t 10 -exclude Glycine_max -generate 
# The final reference sequences can be found in the 353_ref_Fabaceae/353gene directory after downloading

## Explanation of parameters
-o: specifies the output directory
-c: specifiles the taxonomy of species used as reference
-t: specifies the number of threads used to download files
-exclude: excludes species that are not used as reference
-generate: generates a CSV that records the information about the downloaded species
```

2. Run the Easy353, the following is an example of the command syntax:

```shell
easy353.py -1 <input_file1> -2 <input_file2> -r <reference_dir> -o <output_dir> -k1 <filter_kmer> -k2 <assemble_kmer> -t1 <filter_thread> -t2 <assemble_thread>

## The parameters used in the command are explained below:
-1 -2: the input files with paired-end reads, given in FASTQ format. 
-r: the reference directory
-o: the output directory
-k1: K-mer length setting for filtering
-k2: K-mer length setting for assembly
-t1: the threads setting for reads filtering
-t2: the threads setting for reads assembly

# An example
easy353.py -1 Gmax_sim_1.fastq.gz -2 Gmax_sim_2.fastq.gz -r Fabaceae353 -o test_package -k1 31 -k2 41 -t1 1 -t2 4
```

3. Now, you can view the result of Easy353 within the output directory:

The output directory comprises two subdirectories: `filtered_reads` and `target_genes`.

- The `filtered_out` subdirectory contains filtered reads associated with target genes.
- The `assemble_out`subdirectory is the primary result directory of Easy353, storing target gene sequences recovered by the tool.

### Easy353-GUI

Note: There is a detailed guide of Easy353-GUI at 

Please use the Easy353-gui.app or Easy353-gui.exe from https://github.com/plant720/Easy353/releases, the `Easy353-gui.py` is only used for development.

* For Easy353-GUI, you should set the `Paired fq file 1`, `Paired fq file 2`, `Taxonomy` ,and `Output dir` as shown in the image below

![image-20220629211901122](https://cdn.jsdelivr.net/gh/plant720/TyporaPic/img/20220629211901.png)

## User manual

A more complete manual is here: https://github.com/plant720/Easy353/blob/master/docs/manual.pdf

## Contact

If you have any questions, suggestions, or comments about Easy353, feel free to contact the developer at [zzhen0302@163.com](mailto:zzhen0302@163.com).
