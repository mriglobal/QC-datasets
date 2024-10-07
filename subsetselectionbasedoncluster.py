import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import argparse
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import umap
import os
from collections import Counter


def dataframe_to_points(df):
    """Converts 'X' and 'Y' columns of a pandas DataFrame into a NumPy array."""
    if 'X' in df.columns and 'Y' in df.columns:
        points = df[['X', 'Y']].to_numpy()
        return points
    else:
        raise ValueError("DataFrame must contain 'X' and 'Y' columns.")

def main():
    global df, df_cluster, df_taxid_cluster,keep_indices 
    parser = argparse.ArgumentParser(description="Load a TSV file and select cluster representatives.")
    parser.add_argument("-f", "--filepath", type=str, required=True, help="The filepath to the TSV file.")
    parser.add_argument("-n", "--num_representatives", type=int, required=True, help="Number of representatives for each taxid in each cluster.")
    parser.add_argument("-s", "--selection_procedure", type=str, default="centroidbasedsetcover", help="Procedure to select cluster representatives.")
    parser.add_argument("-d", "--data_correction_method", type=str, default="drop", help="method to deal with presumably incorrectly labeled genomes.")
    
    args = parser.parse_args()
    
    try:
        df = pd.read_csv(args.filepath, sep='\t')
        #remove -1 taxids - not useful
        print(len(df))
        df = df[df['hierarchical_clusters'] >= 1]
        print(len(df))
        df_ref = df.copy()
        print("Data loaded successfully. Here's a preview:")
        print(df.head())
        #points = dataframe_to_points(df)

    except Exception as e:
        print(f"An error occurred: {e}")

    if 'selected_rank_taxid' in df:selectlab='selected_rank_taxid'
    else:selectlab='taxid'

    idxtoremove=[]
    for cluster_id, df_cluster in df.groupby('hierarchical_clusters'):
        print(f"Processing hierarchical cluster to quality assure correct taxid labeling: {cluster_id}")
        majority_label = df_cluster[selectlab].mode()[0]
        total_num = len(df_cluster)
        majority_num = (df_cluster[selectlab] == majority_label).sum()
        print(f"majority_label:",majority_label,majority_num,total_num)
        
        for taxid_cluster_id, df_taxid_cluster in df_cluster.groupby('rank'):
            print(f"\tProcessing taxid cluster: {taxid_cluster_id} within rank {cluster_id}")
            print(f"\tUnique taxids:",df_taxid_cluster['taxid'].unique())
            if args.data_correction_method=="drop":
                 if majority_num  > total_num/2: # have a true majority
                     # Identify indices of rows with 'selected_rank_taxid' not in the majority group
                     non_majority_indices = df_taxid_cluster[df_taxid_cluster[selectlab] != majority_label].index
                     # Append these indices to the list of indices to remove
                     idxtoremove.extend(non_majority_indices.tolist())
            if args.data_correction_method=="relabel":
                 if majority_num > total_num / 2:  # have a true majority                
                     # Find indices in the cluster where 'selected_rank_taxid' is not the majority
                     non_majority_indices = df_taxid_cluster[df_taxid_cluster[selectlab] != majority_label].index
                
                     # Relabel 'taxid_clusters' for these indices to the majority taxid cluster value
                     df.loc[non_majority_indices, 'taxid'] = majority_label

    print("number of unique assembly accession entries before filtering {}".format(len(set(df['accession-file']))))
    print("number of unique taxids at rolled up level before filtering {}".format(len(set(df[selectlab]))))
    print(sorted(list(set(df[selectlab]))))
    df.drop(idxtoremove, inplace=True)
    print("number of unique assembly accession entries after filtering {}".format(len(set(df['accession-file']))))
    print("number of unique taxids at rolled up level after filtering {}".format(len(set(df[selectlab]))))
    print(sorted(list(set(df[selectlab]))))

    dfonlytaxidfiltering=df.copy()
  
    idxtoremove=[]
    #print(args.num_representatives)
    for taxid_cluster_id, df_taxid_cluster in df.groupby(selectlab):
        print(f"\tProcessing taxid cluster: {taxid_cluster_id}")
        #print(len(df_taxid_cluster), args.num_representatives)
        if len(df_taxid_cluster) > args.num_representatives: 
            # If the cluster size is larger than the number of representatives, perform subsampling
            sampled_indices = df_taxid_cluster.sample(n=args.num_representatives, random_state=1).index
            # Find indices of rows not sampled and mark them for removal
            unsampled_indices = df_taxid_cluster.index.difference(sampled_indices)
            idxtoremove.extend(unsampled_indices.tolist())

    #print(len(df))
    #print(idxtoremove)
    df.drop(idxtoremove, inplace=True)
    #print(len(df))
    print("number of unique assembly accession entries after subsampling {}".format(len(set(df['accession-file']))))
    print("number of unique taxids at rolled up level after subsampling {}".format(len(set(df[selectlab]))))
    print(sorted(list(set(df[selectlab]))))
 



    #return df,df_cluster,df_taxid_cluster,keep_indices
    df['taxid_clusters'] = [str(taxid) if taxid >= 1 else str(-1) for taxid in list(df['taxid'])]
    dfonlytaxidfiltering.to_csv(args.filepath.replace(".tsv","_onlytaxidfilt_subsampled.tsv"), sep='\t', index=False)
    df.to_csv(args.filepath.replace(".tsv","_taxidandrand_subsampled.tsv"), sep='\t', index=False)
    print("original number of rows:",len(df_ref))
    #print("how many rows got kept:",Counter(list(df["keep"])))
    print("how many rows per selected taxid were kept",Counter(list(df["taxid_clusters"])))
    all_data = px.scatter(df_ref, x='X',y='Y', color='hierarchical_clusters',symbol='rank',hover_data=['description','accession-file',selectlab],width=800,height=800,opacity=0.2)
    #df_subset=df[df['keep']=="Yes"]
    overlay_data = px.scatter(df, x='X',y='Y', color='hierarchical_clusters',symbol='rank',hover_data=['description','accession-file',selectlab],width=800,height=800,opacity=1)
    for trace in overlay_data.data:
        #print(trace)
        all_data.add_trace(trace)

    print(df.head())
    print(df_ref.head())

    all_data.write_html(args.filepath.replace(".tsv","")+"_overlayed_subsample_cluster_umap.html")
    #df.to_csv(args.filepath.replace(".tsv","")+"_cluster_embed_CORE_subsampled.tsv", sep='\t', index=False)

    


#items=[]
if __name__ == "__main__":
    main()


