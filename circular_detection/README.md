## circular_detection

### contents:

*[usage]()
*[description]()
*[options]()
*[bugs]()
*[citation]()

### usage:

```bash
perl circular_detection.pl [-h] [-v] [-f FASTA] [-k KMER] [-c CADRE]
```

### description:

This script detect circular sequence by looking for exact identical k-mer at the two ends on a cadre sequence of the sequences prodvide in fasta file. In order to be able to predict genes spanning the orgin of circular circular, the first 1,000 nucleotides of each circular contigs are dulicated and added at the contig's end. 

### options: 

```bash
-h, --help	produces the help file
-v, --verbose	boolean option to print out warnings during execution. Warnings and errors are redirected to STDERR. Defaults to no verbose (silent mode).
-f, --fasta	sequence fasta file
-k, --kmer	motif length detected to identify circular sequence (default: 10)
-c, --cadre	length of the inspect fragment in 5' for kmer identity finding.
		If 0 screen all the sequence (default: 0)
```

### bugs

* Submit problems or requests here: https://github.com/meb-team/Tools/issues/

### citation

Written by Corentin Hochart (corentin.hochart@uca.fr), UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE). Released under the terms of the GNU General Public License v3. preprocess version v0.1.0.
- 
