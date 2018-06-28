## mga_parser

### usage:

```bash
perl affiliation.pl [-t file] [-a abundance] [-lca taxonomy] [-env] [-longer]
```

### description:

This script produce a lowest common ancestor affiliation (LCA) of each sequence and the abundance of each lineage from a tabular file.

### options: 

```bash
-h, --help	produces the help file
-v, --verbose	boolean option to print out warnings during execution. Warnings and errors are redirected to STDERR. Defaults to no verbose (silent mode).
-f, --force, force the script by ERASED early project. 
-t,	--taxo,  File with complete taxonomy for each contigs.
-d, --delimiter, Delimiter between the sequence name and features id. Must be dot, pipe or underscore (default: pipe)
-lca, Output file with the result of the LCA.
-a, Output file with the abundance of each lineage. 
--env, Boolean option to not take account of 'unclultured viruses' into LCA process. 
--longer, Keep the longest lineage as possible during LCA process.  
```

### bugs

* Submit problems or requests here: https://github.com/meb-team/Tools/issues/

### citation

Written by Corentin Hochart (corentin.hochart@uca.fr), UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE). Released under the terms of the GNU General Public License v3. preprocess version v0.1.0.
