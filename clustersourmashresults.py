import plotly.io as pio
import plotly.express as px
import argparse
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import umap
import os

def load_dist_matrix(file_path):
    # Assuming the similarity matrix is stored in a CSV or similar format
    # Adjust the loading mechanism based on your actual file format (e.g., CSV, TSV, etc.)
    return np.load(file_path)

def cluster_sequences(d_matrix, threshold=1.0):
    # Ensure the distance matrix is symmetrical and diagonal elements are 0
    distance_matrix = squareform(d_matrix, checks=False)

    # Hierarchical clustering
    Z = linkage(distance_matrix,'ward')

    # Form flat clusters
    clusters = fcluster(Z, t=threshold, criterion='distance')

    return(clusters)

def init_umap(d_mat,d_mat_f,tm_file,initview):
    assert os.path.exists(d_mat_f+".labels.txt"), f"File does not exist: {d_mat_f}.labels.txt"
    
    #with open(d_mat_f+".labels.txt") as label_file:
    #labels = [line.strip() for line in label_file.readlines()]

    taxid_table = pd.read_table(tm_file,sep=",",index_col=False)
    print(taxid_table)
    labels=[]
    with open(d_mat_f+".labels.txt", 'r') as file:
        for line in file:
            # Strip newline characters from the end of the line
            labels.append(line.strip())

    taxid_table['accession'] = taxid_table['description'].str.split(' ',expand=True)[0]

    if initview:manifold = umap.UMAP()

    test = manifold.fit_transform(d_mat)

    embedding = pd.DataFrame(test,columns=['X','Y'])


    embedding['accession'] = [l.split('/')[1] for l in labels]
    taxid_table.set_index('assembly', inplace=True,drop=False)
    embedding.set_index('accession', inplace=True,drop=False)
    #print(embedding,taxid_table)
    #print(embedding["accession"])
    #print(taxid_table["assembly"])
    #print(list(taxid_table.index))
    #print(list(embedding['accession'].index))
    embedding['taxid'] = taxid_table.loc[embedding['accession']]['taxid'].astype(str).values
    embedding['rank'] = taxid_table.loc[embedding['accession']]['rank'].values
    embedding['description'] = taxid_table.loc[embedding['accession']]['description'].values
    embedding['parent_taxid'] = taxid_table.loc[embedding['accession']]['parent_taxid'].values
    embedding.rename(columns={'accession': 'accession-file'}, inplace=True)
    print(embedding,taxid_table)
    print(len(embedding),len(taxid_table))
    if initview:fig=px.scatter(embedding, x='X',y='Y', color='taxid',symbol='rank',hover_data=['accession-file','description','parent_taxid'],width=800,height=800)
    else: fig=None
    #fig.show()
    return(embedding,fig)




def main():
    parser = argparse.ArgumentParser(description='Cluster sequences based on a similarity matrix.')
    parser.add_argument('-d', '--dist_matrix', required=True, 
                        help='Path to the file containing the similarity matrix.')
    parser.add_argument('-m', '--taxonomy_metadata', required=True,
                        help='Path to the file containing the similarity matrix.')
    parser.add_argument('-t', '--threshold', type=float, required=True, 
                        help='Similarity threshold for clustering.')
    parser.add_argument('-i', '--inital_view', required=False, action='store_true', 
                        help='flag if you wish to see the umap of the original labels by taxid')
    

    args = parser.parse_args()

    # Load the similarity matrix from the specified file
    dist_matrix = load_dist_matrix(args.dist_matrix)
    print(dist_matrix.shape)

    embed,fig_init=init_umap(dist_matrix,args.dist_matrix,args.taxonomy_metadata,args.inital_view)
    #pio.write_image(fig_init, 'test.png')
    if args.inital_view:fig_init.write_html("initial_clusters_with_taxid.html")
    #fig_init.show()
 
    # Perform clustering with the loaded similarity matrix and the specified threshold
    clusters=cluster_sequences(dist_matrix, args.threshold)

    embed['hierarchical_clusters'] = [str(c) for c in clusters]
    print(embed)
    fig_final = px.scatter(embed, x='X',y='Y', color='hierarchical_clusters',symbol='rank',hover_data=['description','accession-file','taxid'],width=800,height=800)
    fig_final.write_html("taxonomy_rectified_clusters.html")
    embed.to_csv('taxonomic_sequence_cluster.tsv', sep='\t', index=False)





if __name__ == "__main__":
    main()

