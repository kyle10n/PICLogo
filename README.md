# PhyInC Logo: Phylogenetically Independent Contrasts Sequence Logo

PhyInC Logo (pronounced "Fink" for short) is a tool to take Fasta files and Newick tree files to create
a phylogenetically conditioned sequence logo, using Felsensteins _phylogenetic independent
contrast_ (PIC).

## Install from PyPI.org

```
pip3 install phyinc
```
It is always a good idea to use virtual environments, see below.

## Install from a downloaded GitHub repository

For dependencies, start a virtual env as good practice:
```
> python3 -m venv .
> source venv/bin/activate
> pip3 install biopython weblogo matplotlib numpy
```

If you are on mac
```
> brew install ghostscript
```

To install the current package:
```
> pip3 install .
```

## Usage

The basic usage is as follows.
```
> phyinc treefile fastafile
```
Since no outfile is given, the logo is output to `fastafile.fa_logo.pdf`. 
You can decide outputfile and format using the `-o` option:
```
> phyinc -o the_logo.png treefile fastafile
```


### Examples


There is example data in the github repository, and there you can 
run this command:
``` 
> phyinc examples/synthetic_data/ex1_t1.tree examples/synthetic_data/ex1.fa
```
This should create a PDF named "ex1.fa_seqlogo.pdf" in the examples folder.

There are cases when a phylogeny is created on proteins but the logo is
created for domains and this may cause protein accessions be something like 
`ETA_STAAU` but the domain accession is `ETA_STAAU/96-110`. You can then use 
the `--coords` option to ignore the domain coordinates when mapping the domain 
to a tree. 
```
> phyinc --coords examples/PF000672.fa examples/PF000672.treefile
```
However, `phyinc` will report an error if there more than one domain per protein:
```
> phyinc --coords PS00027.fa PS00027.treefile
Error: 'ZFH2_DROME' is a protein appearing twice, probably because you have two 
domains from the same protein in the input. If so, you must submit a tree inferred
on the domain sequences, not on the proteins.
```


# License 

GPLv3


# Authors

* Haolin Guo wrote the basic code as part of his BSc project.
* Kyle Tenn helped make the code into a Python/PyPI package.
* Lars Arvestad oversaw and helped with creating the final Python/PyPI package.

