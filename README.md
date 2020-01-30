# Example: creating new multivec files with row metadata

Set up and activate conda environment:
```
conda env create -f environment.yml
source activate clodius-multivec-env
```

### Notes:

This repo demonstrates the creation of a new multivec file from epilogos chromatin states data<br>
**Goal**: specify row metadata as JSON -> access metadata objects in HiGlass


#### Steps

1. Download and uncompress epilogos file [epilogos/chr1_127epigenomes_15observedStates.txt.gz](https://github.com/Altius/epilogos/blob/master/data/chr1_127epigenomes_15observedStates.txt.gz)
2. Generate row-infos file (containing one JSON object per row)
3. Run clodius `bedfile-to-multivec`  commands [clodius/multivec.rst](https://github.com/higlass/clodius/blob/develop/docs/genomic/multivec.rst#other-data-multivec)

#### Example

```
snakemake
```


### Older notes:

```sh
# without row_info file
clodius convert bedfile-to-multivec \
        chr1_127epigenomes_15observedStates.txt \
        --assembly hg19 \
        --starting-resolution 200
```

To generate row-infos file:

```py
import pandas as pd
df = pd.read_csv('chr1_127epigenomes_15observedStates.txt', sep='\t', header=None)
r = df.columns.values.tolist()[3:]
r = ["State %d" % x for x in r]
s = pd.Series(r)
s.to_frame().to_csv('rowinfo.tsv', sep='\t', index=False, header=False)
```

```
# with row_info file
clodius convert bedfile-to-multivec \
        chr1_127epigenomes_15observedStates.txt \
        --assembly hg19 \
        --starting-resolution 200 \
		--row-infos-filename rownames.tsv \
    	--num-rows 127
```

Output file: `chr1_127epigenomes_15observedStates.multires.mv5`

To verify that row names are included:

```py
import h5py
f = h5py.File('chr1_127epigenomes_15observedStates.multires.mv5', 'r')
f['resolutions']['200'].attrs.keys()
f['resolutions']['200'].attrs['row_infos']
```

Viewconf with horizontal-multivec (with row names) as the center track http://higlass.io/l/?d=e-gJlHZPQ5S3GjCHA86ZMA