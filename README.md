#freqseq

##Intro
`freq_seq` is a small c++ executable that uses the `bamtools` api to record
the frequencies of major and minor 'alleles' from reads mapped in a bam file. 

The program was written a speific task (detecting multiple infectoins from NGS
data), and thus many of the interesting paramters (minor allele frequency cut
off, mapping and base quality cuts) are current hard-coded (though this can be
changed very easily). 

##Build

The build is managed by cmake. This series of commands shuold get you an
exectuable called `freq_seq` :

```
mkdir buid
cd build
cmake ../
make
```

##Usae 

```
./freq_seq [path_to_BAM] [out_prefix]
```

Will produce two files
* `out_prefix_bases.tsv` which contains the minor allele frequency for each
   sample for every base that is polymorphic, or has at least one sample with a
   minor allele at frequency > 0.05
* `out_prefix_summary.tsv` which lists proportion of sites with coverage > 10x
   and a minor allele frequency > 0.2 across the whole genome.


