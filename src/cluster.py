import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage

if __name__ == "__main__":
    print(snakemake.input[0])

    df = pd.read_csv(snakemake.input[0], sep='\t', header=None)
    colnames = df.columns.values.tolist()
    print(colnames)