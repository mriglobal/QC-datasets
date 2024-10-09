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
```
python QC_dataset.py -i input_data/Francisella_ncbi.zip -bm median -bp .1 -n 3 -d --gc_lower .3 --gc_upper .35 -se plasmid -on f_dataset_forqc_workflow -of directory -sd
```
filter inital set of genomes for quality. remove records that have too many Ns, have a weird %GC, or a label that would indicate bad quality or limited use for comparing different genomes like 'partial', 'hypothetical', 'plasmid', etc

```
usage: QC Check - provide datasets zip, a multi fasta, or a directory of fasta files. Perform QC as desired. Output one of the three formats
       [-h] -i INPUT [-n N_FILTER] [-d] [--gc_lower GC_LOWER] [--gc_upper GC_UPPER] [-bm BIN_METHOD] [-bp BIN_PERCENTILE]
       [-si [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]]] [-se [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]]] -on OUTPUT_NAME [-of {fasta,directory,zip}] [-sd] [-r]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        NCBI Datasets zip file, multifasta file, or directory of fasta files to read in
  -n N_FILTER, --n_filter N_FILTER
                        input criteria for filtering by N characters in sequence. If a float, filters by that percent of Ns in sequence. If an integer, filters any
                        sequences with that many consecutive Ns. If "any", filter records containing any number of N
  -d, --degen_filter    flag if all records containing any degen characters should be filtered
  --gc_lower GC_LOWER   lower bound of GC content to be kept in sequences
  --gc_upper GC_UPPER   upper bound of GC content to be kept in sequences
  -bm BIN_METHOD, --bin_method BIN_METHOD
                        input criteria for filtering by N characters in sequence. Either "median" if binning should be done by median sequence length, or specify an
                        integer length
  -bp BIN_PERCENTILE, --bin_percentile BIN_PERCENTILE
                        float value for percent difference in length to cut off binning (+- 10 percent by default)
  -si [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]], --strings_to_include [STRINGS_TO_INCLUDE [STRINGS_TO_INCLUDE ...]]
                        An optional space-separated list of strings to filter inclusively
  -se [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]], --strings_to_exclude [STRINGS_TO_EXCLUDE [STRINGS_TO_EXCLUDE ...]]
                        An optional space-separated list of strings to filter exclusively
  -on OUTPUT_NAME, --output_name OUTPUT_NAME
                        What name should the output have
  -of {fasta,directory,zip}, --output_format {fasta,directory,zip}
                        Output format: fasta, directory, or zip.
  -sd, --save_drops     true or false to write out all filtered sequences
  -r, --refseq_only     keep only refseq type headers in output
```

#Stage 2 - setting up filter folder for further analysis

only need to do this if -of is different than directory, and will need to process the output so that it looks like the output of directory

#Stage 3A - use sourmash to produce a distance matrix (.npz file), that can be used with umap to cluster records based on sequence properties

after making sure sourmash is working correctly, get initial distance matrix to be used later in clustering

```
cd f_dataset_forqc_workflow/
for fil in *fna;do sourmash sketch dna -p k=31,scaled=100,noabund $fil ; done
cd ..
sourmash compare --processes 20 --distance-matrix --ksize 31 --dna --scaled 100 f_dataset_forqc_workflow/*.sig -o f-dataset-comparison-matrix.npz
cd $fol
```
see sourmash documentation for more details if you wish the understadn all the options in makind the distacen matrix

#Stage 3B - taxonomy roll up for records, to allow for labeling by taxid, rank, etc.

here is where the taxonomy metadata is labeled onto the seqs

```
time python gettaxonomyfororgs.py -d f_dataset_forqc_workflow/
```
```
usage: gettaxonomyfororgs.py [-h] -d DIRECTORY

Process FASTA files in a folder of folders.

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        The path to the directory of directories containing FASTA files.
```


#Stage 4 -  after all preprocessing is done, use the inforamtion to generate a umap representation of clustering

Doing an init run of this script can also initial view the records clustered by initial taxid each accession is assigned to, 
always will then make a umap of the records, using the distances provided by sourmash and produced by cutting clusters that differ more than the treshold -t

Note: for now always do a init run with the -i flag

```
python clustersourmashresults.py -d f_dataset_forqc_workflow/ncbi_dataset/data/f-dataset-comparison-matrix.npz -m f_dataset_forqc_workflow_taxonomy.csv  -t 1.0 -i
```
```
usage: clustersourmashresults.py [-h] -d DIST_MATRIX -m TAXONOMY_METADATA -t THRESHOLD [-i]

Cluster sequences based on a similarity matrix.

optional arguments:
  -h, --help            show this help message and exit
  -d DIST_MATRIX, --dist_matrix DIST_MATRIX
                        Path to the file containing the similarity matrix.
  -m TAXONOMY_METADATA, --taxonomy_metadata TAXONOMY_METADATA
                        Path to the file containing the similarity matrix.
  -t THRESHOLD, --threshold THRESHOLD
                        Similarity threshold for clustering.
  -i, --inital_view     flag if you wish to see the umap of the original labels by taxid
```
#Stage 5 - run subsampling to maintain as much variability in the sequences while minimizing total sequence amount retained

```
python  data_utils/workflowforatds/subsetselectionbasedoncluster.py -f taxonomic_sequence_cluster.tsv -n 20
```
```
usage: subsetselectionbasedoncluster.py [-h] -f FILEPATH -n NUM_REPRESENTATIVES [-s SELECTION_PROCEDURE] [-d DATA_CORRECTION_METHOD]

Load a TSV file and select cluster representatives.

optional arguments:
  -h, --help            show this help message and exit
  -f FILEPATH, --filepath FILEPATH
                        The filepath to the TSV file.
  -n NUM_REPRESENTATIVES, --num_representatives NUM_REPRESENTATIVES
                        Number of representatives for each taxid in each cluster.
  -s SELECTION_PROCEDURE, --selection_procedure SELECTION_PROCEDURE
                        Procedure to select cluster representatives.
  -d DATA_CORRECTION_METHOD, --data_correction_method DATA_CORRECTION_METHOD
                        method to deal with presumably incorrectly labeled genomes.

```

