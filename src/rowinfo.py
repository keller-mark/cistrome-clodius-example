import pandas as pd
import json
import random

random.seed(10)

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0], sep='\t', header=None)
    state_cols = map(str, df.columns.values.tolist()[3:])

    with open(snakemake.input[1], 'r') as f:
        tree_matrix = json.load(f)
    
    with open(snakemake.output[0], 'w') as f:
        for state_col in state_cols:
            metadata = {
                "state": state_col,
                "h": tree_matrix[state_col],
                "r1": random.randint(0, 10),
                "r2": random.randint(0, 100)
            }
            metadata_json = json.dumps(metadata)
            f.write(f"{metadata_json}\n")