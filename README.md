# Example: creating new multivec files with row metadata

Set up and activate conda environment:
```
conda env create -f environment.yml
source activate clodius-multivec-env
```

Run the example:

```
snakemake
```

### Notes:

This repo demonstrates the creation of a new multivec file from epilogos chromatin states data<br>
**Goal**: specify row metadata as JSON -> access metadata objects in HiGlass

Output file: `chr1_127epigenomes_15observedStates.multires.mv5` in `data/processed/`

To verify that row names are included:

```py
import h5py
f = h5py.File('chr1_127epigenomes_15observedStates.multires.mv5', 'r')
f['resolutions']['200'].attrs.keys()
f['resolutions']['200'].attrs['row_infos']
```

Viewconf with horizontal-multivec (with row names) as the center track http://higlass.io/l/?d=e-gJlHZPQ5S3GjCHA86ZMA