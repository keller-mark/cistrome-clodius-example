import json
from uuid import uuid4

# Encode a tree of nested JSON objects as a "matrix"
def tree_to_matrix(tree):

    def is_leaf(node):
        return ("children" not in node.keys() or len(node["children"]) == 0)

    used_uuids = []

    def convert(node, prev):
        curr = prev.copy()

        # Base case:
        if is_leaf(node):
            curr.append(node["name"])
            return curr
        
        # General case:
        level_id = str(uuid4())
        assert(level_id not in used_uuids)
        used_uuids.append(level_id)
        curr.append(level_id)

        # Recursively run `convert()` on the node's children:
        children = map(lambda c: convert(c, curr), node["children"])
        result = []
        for child in children:
            if isinstance(child, list) and isinstance(child[0], list):
                result = result + child
            else:
                result.append(child)
        
        return result
    
    matrix = convert(tree, [])

    # `matrix` is a 2-dimensional list, where each inner array represents the position of a leaf node.

    # We want to convert the uuid strings to unique integers to reduce the size of the data.
    # We also want to map each leaf to its array, so we create `matrix_dict`
    #   where key is leaf name, value is the array associated with the leaf.
    matrix_dict = dict()
    uuid_to_int = dict(zip(used_uuids, map(lambda i: f"i-{i}", range(len(used_uuids)))))
    for i in range(len(matrix)):
        for j in range(len(matrix[i]) - 1):
            matrix[i][j] = uuid_to_int[matrix[i][j]]
        
        matrix_dict[matrix[i][len(matrix[i]) - 1]] = matrix[i]

    return matrix_dict
    

if __name__ == "__main__":
    with open(snakemake.input[0], 'r') as f:
        tree_dict = json.load(f)

    tree_matrix_dict = tree_to_matrix(tree_dict)

    with open(snakemake.output[0], 'w') as f:
        json.dump(tree_matrix_dict, f)
