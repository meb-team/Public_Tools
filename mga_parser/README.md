## mga_parser

### usage:

```bash
 perl mga_parser.pl [-i file.fasta] [-aa file.faa] [-nt file.fna] [-tsv file.tsv] [-gff file.gff]
```

### description:

This script automatically predict gene/protein from a fasta file with [MetaGeneAnnotator](http://metagene.nig.ac.jp/) predictor tool (Hideki Noguchi, Takeaki Taniguchi and Takehiko Itoh, DNA Research, 2008). 

### options: 

```
-h, --help	produces the help file
-v, --verbose	boolean option to print out warnings during execution. Warnings and errors are redirected to STDERR. Defaults to no verbose (silent mode).
-f, --force, force the script by ERASED early project. 
-i, --input, nucleotid fasta file as input.
--procedure [m|s]	
	Multi or Single species option in MetageneAnnotator. 
	Meta or single parameter for -p parameter in prodigal (default 'm').
--mga, call existing mga file for nucleotid fasta file (-i).
--nt, output predicted gene in fasta format.
--aa, output predicted protein in fasta format.
--gff, output predicted gene in gff format.
--tsv, output predicted protein and gene in tabular format. 
```

### bugs

* Submit problems or requests here: https://github.com/meb-team/Tools/issues/

### citation

Written by Corentin Hochart (corentin.hochart@uca.fr), UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE). Released under the terms of the GNU General Public License v3. preprocess version v0.1.0.
