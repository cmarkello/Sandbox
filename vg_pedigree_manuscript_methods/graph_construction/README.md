## Running the scripts

First, download the requisite vcf files and generate the segmental duplicate-free vcfs.

```
./preprocess_1000GPlons_GRCh38_graph.sh
```

Next, run the construction workflow for building the graph and indexes for the primary graph reference.

```
./make_primary_GRCh38_graph.sh
```

Finally, run the construction workflow for building the graph and indexes for the non-segmental-duplicate 1000GP-based graph reference.

```
./make_1000GPlons_GRCh38_graph.sh
```

## Output files

