import argparse
import pandas as pd
import os
from Bio import SeqIO
from MRItaxonomy import accession2taxid,taxid
def_ba_taxid=accession2taxid.get_taxid('NC_007530.2')
def_ba_rank=taxid.getrank(def_ba_taxid)
import time
import argparse
import pandas as pd
import os
from Bio import SeqIO

def process_fasta_file(filepath):
    """
    Reads headers from a FASTA file, splits each header into two parts,
    and saves the first part while keeping the second part as a list.

    Parameters:
    - filepath: Path to the FASTA file.

    Returns:
    - A pandas DataFrame with 'assembly' and 'header_part_2' columns.
    """
    data = []
    for record in SeqIO.parse(filepath, "fasta"):
        #print(record.description)
        header_parts = record.description.split(" ", 1)  # Assuming space separates two parts in the header
        #print(header_parts)
        if len(header_parts) == 2:
            data.append({'accession': header_parts[0], 'assembly':filepath.split("/")[-1], 'description': header_parts[1]})
        else:
            print(f"Header format unexpected in file {filepath}: {record.id}")
    
    return pd.DataFrame(data)

def process_directory(directory):
    """
    Walks through a directory of directories to process each FASTA file.

    Parameters:
    - directory: Path to the directory containing subdirectories with FASTA files.
    """
    df_list=[]
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fna"):
                filepath = os.path.join(root, file)
                df_list.append(process_fasta_file(filepath))
                if not df_list[-1].empty:
                    print(f"metadata saved to dataframe:",filepath)
                else:
                    print("no metadata harvested from:",filepath)
    combined_df = pd.concat(df_list, ignore_index=True)
    #combined_df['taxid']=[accession2taxid.get_taxid(item) for item in list(combined_df['assembly'])]
    #combined_df['rank']=[taxid.getrank(item) for item in list(combined_df['taxid'])]
    combined_df.insert(len(combined_df.columns)-1, 'taxid', [addprintout(accession2taxid.get_taxid(item)) for item in list(combined_df['accession'])], allow_duplicates=True)
    combined_df.insert(len(combined_df.columns)-1, 'rank', [addprintout(taxid.getrank(item),"taxid.getrank") for item in list(combined_df['taxid'])], allow_duplicates=True)
    combined_df['parent_taxid']=[addprintout(taxid.getparent(item),"taxid.getparent") for item in list(combined_df['taxid'])]
    combined_df['parent_rank']=[addprintout(taxid.getrank(item),"taxid.getrank") for item in list(combined_df['parent_taxid'])]
    #combined_df['accession']=[taxid.getrank(item) for item in list(combined_df['parent_taxid'])]

    folder=directory.split('/')[0]
    output_filename = os.path.join(os.getcwd(),folder+"_taxonomy.csv")
    print(combined_df,output_filename)
    combined_df.to_csv(output_filename, index=False)

def addprintout(val,string="accession2taxid.get_taxid"):
    print(string, time.time())
    return(val)


def main():
    parser = argparse.ArgumentParser(description="Process FASTA files in a folder of folders.")
    parser.add_argument("-d", "--directory", type=str, required=True, help="The path to the directory of directories containing FASTA files.")

    args = parser.parse_args()

    process_directory(args.directory)

if __name__ == "__main__":
    main()

