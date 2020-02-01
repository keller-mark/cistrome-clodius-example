# Example: creating new multivec files with row metadata

Set up and activate conda environment:

```sh
conda env create -f environment.yml
source activate clodius-multivec-env
```

Run the example:

```sh
snakemake
```

Load processed multivec file into [higlass server](https://github.com/higlass/higlass-server):

```sh
# source deactivate
# cd path/to/higlass-server
# workon higlass-server
python manage.py ingest_tileset \
    --uid cistrome-demo \
    --filename path/to/data/processed/chr1_127epigenomes_15observedStates.multires.mv5 \
    --filetype multivec \
    --datatype matrix
# python manage.py runserver 9000
# deactivate
```

Use the following track definition:

```json
{
    "type": "horizontal-multivec",
    "uid": "cistrome-track",
    "tilesetUid": "cistrome-demo",
    "server": "http://localhost:9000/api/v1",
    "options": {
        "labelPosition": "hidden",
        "labelColor": "black",
        "labelTextOpacity": 0.4,
        "valueScaling": "linear",
        "trackBorderWidth": 0,
        "trackBorderColor": "black",
        "heatmapValueScaling": "log",
        "name": "chr1_127epigenomes_15observedStates.multires.mv5",
        "labelLeftMargin": 0,
        "labelRightMargin": 0,
        "labelTopMargin": 0,
        "labelBottomMargin": 0,
        "labelShowResolution": true,
        "minHeight": 100
    },
    "width": 1607,
    "height": 482,
    "resolutions": [
        16384000,
        8192000,
        4096000,
        2048000,
        1024000,
        512000,
        256000,
        128000,
        64000,
        32000,
        16000,
        8000,
        4000,
        2000,
        1000
    ]
}
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