# cd-CAP

### System Requirements
- make (version 3.81 or higher)
- g++ (GCC version 4.1.2 or higher)
- IBM ILOG CPLEX Optimization Studio

### Compiling cd-CAP
In the `Makefile`, set `CPLEXROOT` to the path of your root CPLEX folder.

The Makefile is set-up for GCC 6.2. If you are using GCC version 4.x, add ` std=gnu++0x` flag to `CCC` in the Makefile.

Simply run `make` command in the root cd-CAP folder. It will create the executables.

### Running `mcsc`
**Usage:**
```sh
./mcsc -n [network] -l [alteration profiles] -c [chromosome information; optional] -r [min number of colours in subnetwork] -x [exclude genes; optional] -s [maximum subnetwork size] -t [minimum subgraph recurrence]  -k [number of subnetworks] -e [error; optional] -f [outputFolder] -d [threads] -t [time limit in seconds]
```

| Parameters | Description |
| ------ | ------ |
| `-n` | network file | 
| `-l` | alterations file |
| `-c` | (optional) gene-to-chromosome map |
| `-x` | Excluded genes |
| `-f` | output folder name |
| `-r` | minimum number of colors in each subnetwork |
| `-s` | maximum subnetwork size |
| `-t` | minimum sample recurrence |
| `-k` | number of resulting subnetworks |
| `-e` | (optional) allowed extension error rate |
| `-d` | number of threads used |
| `-h` | time limit in seconds |



`-n` : &nbsp;&nbsp; This parameter represents an edge collection file where each row represents an edge in form of two node names, separated by whitespace. All edges are treated as undirected. There is no header row. e.g.
```sh
	A1BG    CRISP3
	A1CF    APOBEC1
	A2M     ABCA1
	...
```

`-l` : &nbsp;&nbsp; This parameter represents a file containing information about alterations in all the input samples, in form of "SampleID Gene AltType" rows. Currently, up to 64 different alteration types are supported (the third column). There is no header row. e.g.
```sh
	T294    CCNL2   SNV
	T294    PTCHD2  SNV
	T294    COL16A1 SNV
	...
```

`-c` : &nbsp;&nbsp; This `optional` parameter is a gene-chromosome map, allowing for more information in the output. 
```sh
	Gene    Chromosome      KaryotypeBand
	ADAM30  1       p12
	HAO2    1       p12
	HMGCS2  1       p12
	...
```

`-x` : &nbsp;&nbsp; This `optional` parameter represents a file containing a list of genes whose colors should be removed after reading the input.
```sh
	Gene1
	Gene2
	Gene3
	...
```

`-f` : &nbsp;&nbsp; This parameter contains the name of the output folder for the run of the program, in which all the output files will be stored. This folder name will have values of parameters below appended to it.

`-r` : &nbsp;&nbsp; This integer parameter controls the minimum required number of colors among the nodes of each resulting subnetworks, i.e. how "colorful" a subnetwork must be. Keep the value set to 1 for default configuration.

`-s` : &nbsp;&nbsp; This integer parameter controls the maximum subnetwork size. For the first time running the program on a new dataset, 10 could be a reasonable value.

`-t` : &nbsp;&nbsp; This integer parameter controls the minimum required sample recurrence of each resulting subnetwork.

`-k` : &nbsp;&nbsp; This integer parameter controls the number of subnetworks that we wish to detect.

`-e` : &nbsp;&nbsp; This `optional` floating type parameter controls the maximum allowed error rate when extending subnetworks before the optimization. If not specified, it defaults to 0.

`-t` : &nbsp;&nbsp; This integer parameter controls the minimum required sample recurrence of each resulting subnetwork.

`-d` : &nbsp;&nbsp; This integer parameter specifies the number of threads used for the optimization.

`-h` : &nbsp;&nbsp; This integer parameter specifies the number of seconds that the optimization step is allowed to take before returning a solution.

#### Example
`./mcsc -n ../data/STRING10_HiConf_PPI.edges -l ../data/alteration_status_COAD_20171108.tsv -c ../data/string10_node_chromosome_map.tsv -r 1 -s 10 -t 138 -k 100 -e 0 -d 32 -h 36000 -f TCGA_COAD`
