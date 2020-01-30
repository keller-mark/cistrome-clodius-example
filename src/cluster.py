import pandas as pd
import numpy as np
import scipy.cluster
from functools import reduce
import json

# Create a nested dictionary from the ClusterNode's returned by SciPy
def add_node(node, parent):
    # First create the new node and append it to its parent's children
    new_node = dict( node_id=node.id, children=[] )
    parent["children"].append( new_node )
    # Recursively add the current node's children
    if node.left: add_node( node.left, new_node )
    if node.right: add_node( node.right, new_node )

def cluster(df):
    state_cols = df.columns.values.tolist()[3:]
    state_df = df[state_cols].transpose()
    observation_vectors = state_df.values
    observation_labels = list(state_df.index.values)

    Z = scipy.cluster.hierarchy.linkage(observation_vectors, method='ward')
    T = scipy.cluster.hierarchy.to_tree(Z)

    # Create dictionary for labeling nodes by their IDs
    id_to_label = dict(zip(range(len(observation_labels)), observation_labels))

    # Initialize nested dictionary, then recursively iterate through tree
    tree_dict = dict( children=[], name="root" )
    add_node(T, tree_dict)

    # Label each node with the names of each leaf in its subtree
    def label_tree( n ):
        # If the node is a leaf, then we have its name
        if len(n["children"]) == 0:
            leaf_names = [ id_to_label[n["node_id"]] ]
        # If not, flatten all the leaves in the node's subtree
        else:
            leaf_names = reduce(lambda ls, c: ls + label_tree(c), n["children"], [])
        # Delete the node id since we don't need it anymore and it makes for cleaner JSON
        del n["node_id"]
        # Labeling convention: "-"-separated leaf names
        n["name"] = name = "-".join(sorted(map(str, leaf_names)))
        return leaf_names

    label_tree(tree_dict["children"][0])
    return tree_dict


if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0], sep='\t', header=None)

    tree_dict = cluster(df)

    with open(snakemake.output[0], 'w') as f:
        json.dump(tree_dict, f)
