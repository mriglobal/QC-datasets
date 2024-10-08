# QC-datasets

# CasCADE

Evaluate candidate gRNA for exclusivity

### conda environment setup:
```
conda env create -f CasCADE.yml
conda activate CasCADE
pip install -r requirements.txt
```
*Note that on some systems conda may not load library locations correctly for files like /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29', in that case you simply run ```export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"``` before using the environment to fix the issue

# Cheat-sheet for CasCADE design parameters for given target

| **Design Characteristic** | **dsDNA** | **ssDNA** | **dsRNA** | **ss+RNA** | **ss-RNA** |
| --------------------- | ----- | ----- | ----- | ------ | ------ |
| **--single_strand**       | No    | forward | No  | reverse | forward | 
| **Cas Protien**           | 12a   | 12a   | 13d   | 13d    | 13d    |
| **PAM**                   | TTTV, 5' | TTTV, 5' | No | No  | No     |
| **Scaffold Side**         | 5'    | 5'    | 5'    | 5'     | 5'     |


# Scaffolds
**cas13d scaffold** CAAGUAAACCCCUACCAACUGGUCGGGGUUUGAAAC
**cas12a scaffold** UAAUUUCUACUAAGUGUAGAU

# Design candidate gRNAs

### search for conserved kmers in target group
```
python design/grna_inner_set.py -k 20 --seqs dataset.zip --pseq pamseq --pside 5prime/3prime -t threads -o outprefix
```
* -k is the size of conserved kmer to search for (default: 18)
* --seqs is a directory containing sequence files or an ncbi datasets zip file containing sequence files
* --pseq is the pam sequence i.e. TTTV for Cas12a, Note: TTTV, 5' above includes both pseq and pside arguments
* --pside specifies if the pam occurs at the 5 or 3 prime side of the gRNA
* -t is the number of threads to be used (NOTE: may not work if more threads are specified than records provided. Still under investigation.)
* -o is the outprefix (also impacts naming within output file, ideally is descriptive of target group for record keeping)
* -r optional reference genome (must be specified without output prefix)
* --single_strand [forward, reverse] optionally should be passed forward or reverse depending on if targets should only be found on the forward or reverse strand (i.e. RNA target). Default is to consider both strands
* -c is an optional argument to divide data into chunks before counting, and is useful in memory-limited scenarios
* --pstrand currently does nothing, but in the future will support specifying if the PAM occurs on a particular strand
* -p is the percent variance in reference length for replicon binning, and will be typically not passed and left at the default of .1
* --split_records_size is the size to chunk larger records to in megabasepairs. When on eukaryotic organisms, try specifying --split_records_size 1
* outputs a .tsv of the candidate designs and a .fasta of gRNA_kmers

### OPTIONALLY: find k-mer hypergraph that covers entire sequence set
```
python design/grna_hypergraph_complete.py -k 20 --seqs sequences_to_consider
```
