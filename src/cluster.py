import pandas as pd
import numpy as np
import scipy.cluster
from functools import reduce
import json


def linkage_matrix_to_tree(Z, labels):
    # Reference: https://stackoverflow.com/a/20060322
    
    T = scipy.cluster.hierarchy.to_tree(Z)
    
    # Create dictionary for labeling nodes by their indices
    id_to_label = dict(zip(range(len(labels)), labels))

    # Create a nested dictionary from the ClusterNodes returned by SciPy
    def add_node(node, parent):
        # First create the new node and append it to its parent's children
        name = node.id if node.left or node.right else id_to_label[node.id]
        new_node = dict( name=str(name) )
        if "children" in parent.keys():
            parent["children"].append( new_node )
        else:
            parent["children"] = [ new_node ]
        # Recursively add the current node's children
        if node.left: add_node( node.left, new_node )
        if node.right: add_node( node.right, new_node )

    # Initialize dictionary, then recursively iterate through tree
    tree_dict = dict( name="root", children=[] )
    add_node(T, tree_dict)
    return tree_dict


def cluster(df):
    state_cols = df.columns.values.tolist()[3:]
    state_df = df[state_cols].transpose()
    observation_vectors = state_df.values
    observation_labels = list(state_df.index.values)

    Z = scipy.cluster.hierarchy.linkage(observation_vectors, method='ward')
    return Z, observation_labels
    

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0], sep='\t', header=None)
    num_rows = int(snakemake.wildcards["num_rows"])
    df = df[df.columns.values.tolist()[:3+num_rows]]
    Z, labels = cluster(df)
    tree_dict = linkage_matrix_to_tree(Z, labels)

    with open(snakemake.output[0], 'w') as f:
        json.dump(tree_dict, f)
