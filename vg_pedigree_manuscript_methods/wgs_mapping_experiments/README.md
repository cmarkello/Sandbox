## Running the scripts

First, run the bwa mem mapper on the HG002 and HG005 trio sample set.

```
./bwamem_map.sh
```

Next, run the Illumina DRAGEN platform mapper on the HG002 and HG005 trio sample set.

```
./dragen_map.sh
```

Following that, run the VG Giraffe alignment on the HG002 and HG005 trio sample sets against the 1000GP graph reference 

```
./giraffe_map.sh
```

Finally, run the VG Pedigree alignment workflow on the HG002 and HG005 trio sample sets.
Once, each, using the default deeptrio and deepvariant models.
And another, each using the trained deeptrio and deepvariant models.

```
./giraffe_pedigree_map.sh HG002 default
./giraffe_pedigree_map.sh HG005 default
./giraffe_pedigree_map.sh HG002 trained
./giraffe_pedigree_map.sh HG005 trained
```

## Output files

