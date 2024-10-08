# QC-datasets
Workflow to polish ncbi datasets .zip (and potentially and multifasta or folder of fastas) of undesirable genomes and even subsample for tractibility  

### conda environment setup:
```
conda create -n qc_tool_dev_env_streamlit python=3.8 plotly pandas scipy umap-learn marisa-trie streamlit biopython sourmash -c conda-forge -c bioconda
conda activate qc_tool_dev_env_streamlit
pip install MRItaxonomy
# then intialize the taxonomy database (only need to do once) by
python
from MRItaxonomy import accession2taxid,taxid
def_ba_taxid=accession2taxid.get_taxid('NC_007530.2')
def_ba_rank=taxid.getrank(def_ba_taxid)
```
*Note that on some systems conda may not load library locations correctly for files like /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29', in that case you simply run ```export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"``` before using the environment to fix the issue

# Cheat-sheet for 


# Sample work flow for use

#Stage 1 - basic genome QC filtering
#example command:
```
python QC_dataset.py -i input_data/Francisella_ncbi.zip -bm median -bp .1 -n 3 -d --gc_lower .3 --gc_upper .35 -se plasmid -on f_dataset_forqc_workflow -of zip -sd
```
#filter inital set of genomes for quality. remove records that have too many Ns, have a weird %GC, or a label that would indicate bad quality or limited use for comparing different genomes like 'partial', 'hypothetical', 'plasmid', etc

#usage: QC Check - provide datasets zip, a multi fasta, or a directory of fasta files. Perform QC as desired. Output one of the three formats
#       [-h] -i INPUT [-n N_FILTER] [-d] [--gc_lower GC_LOWER] [--gc_upper GC_UPPER] [-bm BIN_METHOD] [-bp BIN_PERCENTILE]
#       [-si [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]]] [-se [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]]] -on OUTPUT_NAME [-of {fasta,directory,zip}] [-sd] [-r]

#optional arguments:
#  -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        NCBI Datasets zip file, multifasta file, or directory of fasta files to read in
#  -n N_FILTER, --n_filter N_FILTER
#                        input criteria for filtering by N characters in sequence. If a float, filters by that percent of Ns in sequence. If an integer, filters any
#                        sequences with that many consecutive Ns. If "any", filter records containing any number of N
#  -d, --degen_filter    flag if all records containing any degen characters should be filtered
#  --gc_lower GC_LOWER   lower bound of GC content to be kept in sequences
#  --gc_upper GC_UPPER   upper bound of GC content to be kept in sequences
#  -bm BIN_METHOD, --bin_method BIN_METHOD
#                        input criteria for filtering by N characters in sequence. Either "median" if binning should be done by median sequence length, or specify an
#                        integer length
#  -bp BIN_PERCENTILE, --bin_percentile BIN_PERCENTILE
#                        float value for percent difference in length to cut off binning (+- 10 percent by default)
#  -si [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]], --strings_to_include [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]]
#                        An optional space-separated list of strings to filter inclusively
#  -se [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]], --strings_to_exclude [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]]
#                        An optional space-separated list of strings to filter exclusively
#  -on OUTPUT_NAME, --output_name OUTPUT_NAME
#                        What name should the output have
#  -of {fasta,directory,zip}, --output_format {fasta,directory,zip}
#                        Output format: fasta, directory, or zip.
#  -sd, --save_drops     true or false to write out all filtered sequences
#  -r, --refseq_only     keep only refseq type headers in output

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
