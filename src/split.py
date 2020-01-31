import json
from uuid import uuid4

# Encode a tree of nested JSON objects as a "matrix"
def tree_to_matrix(tree):

    def is_leaf(node):
        return ("children" not in node.keys() or len(node["children"]) == 0)

    used_level_ids = []

    def convert(node, prev):
        curr = prev.copy()

        if is_leaf(node):
            curr.append(node["name"])
            return curr
        
        level_id = str(uuid4())
        assert(level_id not in used_level_ids)
        used_level_ids.append(level_id)

        curr.append(level_id)
        children = map(lambda c: convert(c, curr), node["children"])
        result = []
        for child in children:
            if isinstance(child, list) and isinstance(child[0], list):
                result = result + child
            else:
                result.append(child)
        
        return result
    
    matrix = convert(tree, [])

    matrix_dict = dict()
    uuid_to_int = dict(zip(used_level_ids, map(lambda i: f"i-{i}", range(len(used_level_ids)))))
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
