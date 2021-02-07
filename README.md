## Overview ##
Accel-align is a fast alignment tool implemented in C++ programming language.

## Get started ##

### Installation ###
* clone it from code repository
* make it: `make -j`

### Build index ###
It's mandatary to build the index before alignment. 
Options:
```
-l INT the length of k-mers
-m use low memory
```

Example:
```
path-to-accel-align/index -m -l 32 path-to-ref/ref.fna
```
It will generate the index aside the reference genome as `path-to-ref/ref.fna.hash`.

### Align ###
When the alignment trigged, the index will be loaded in memory automaticly.

Options:
```
   -t INT number of cpu threads to use[1].
   -l INT length of seed [32].
   -o name of output file to use.
   -x alignment-free.
   -w use WFA for extension. It's using KSW by default.
   -p the maximum distance allowed between the paired-end reads[1000].
    maximum read length and read name length supported are 512.
```

### Pair-end aligment ### 

``` 
path-to-accel-align/accalign-cpu options ref.fna read1.fastq read2.fastq
```

Example:
``` 
path-to-accel-align/accalign-cpu -l 32 -t 4 -o output-path/out.sam \
path-to-ref/ref.fna input-path/read1.fastq input-path/read2.fastq
``` 

### Single-end aligment ### 

``` 
path-to-accel-align/accalign-cpu options ref.fna read.fastq
```

Example:
``` 
path-to-accel-align/accalign-cpu -l 32 -t 4 -o output-path/out.sam \
path-to-ref/ref.fna input-path/read.fastq
``` 

### Aligment-free mode ### 
Accel-Align does base-to-base align by default. However, Accel-Align supports alignment-free mapping mode where the position is reported without the CIGAR string.
 `-x` option will enable the alignment-free mode.

Example:
``` 
path-to-accel-align/accalign-cpu -l 32 -t 4 -x -o output-path/out.sam \
path-to-ref/ref.fna input-path/read.fastq
``` 
  