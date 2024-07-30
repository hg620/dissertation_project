#Python code used for choosing 500 representative sequences out of 4000
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

def read_distance_matrix(filename):
    matrix = []
    with open(filename, 'r') as file:
        for line in file:
            row = line.strip().split()
            if len(row) > 0:  
                row = [float(x) if x != '-' else np.nan for x in row]
                matrix.append(row)
    return np.array(matrix)

def perform_clustering(distance_matrix, num_clusters=500):
    nan_mask = np.isnan(distance_matrix)
    max_value = np.nanmax(distance_matrix)
    distance_matrix[nan_mask] = max_value
    
    condensed_matrix = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_matrix, method='average')
    cluster_labels = fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    return cluster_labels

def select_representatives(cluster_labels, distance_matrix, sample_names):
    clusters = {}
    for idx, cluster_id in enumerate(cluster_labels):
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(idx)

    representatives = []
    for cluster_id, indices in clusters.items():
        if len(indices) == 1:
            representatives.append(indices[0])
        else:
            submatrix = distance_matrix[np.ix_(indices, indices)]
            representative_idx = np.argmin(np.sum(submatrix, axis=0))
            representatives.append(indices[representative_idx])

    representative_names = [sample_names[idx] for idx in representatives]
    return representative_names

distance_matrix = read_distance_matrix('filtered_output.txt')

with open('sample_names.txt', 'r') as namefile:
    sample_names = [line.strip() for line in namefile.readlines()]

cluster_labels = perform_clustering(distance_matrix, num_clusters=500)
representatives = select_representatives(cluster_labels, distance_matrix, sample_names)

print("Cluster Representatives:", representatives)

with open('representatives.txt', 'w') as file:
    for rep in representatives:
        file.write(f"{rep}\n")
