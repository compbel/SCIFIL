# SCIFIL
Single Cell Inference of FItness Landscape

We propose a computational method for in vivo inference of clonal
selection and estimate of fitness landscapes of heterogeneous cancer cell populations from
single cell sequencing data.

It takes single cell data, mutation tree and estimates finesses of all mutations.

## Parameters
- ``n`` - number of haplotypes
- ``m`` - number of mutation (not counting repeated mutation)
- ``gv_file``  - path to file with mutation tree in GraphcViz format
- ``names_file`` (_optional_) - name of file with mutation names
- ``output`` (_optional_) - path for output file. Default "out.txt"
- ``method`` (_optional_) - method to use for fitness calculation(default "heuristic"):
    - ``heuristic`` to use heuristic
    - ``brute_force`` to use brute force in finding order of mutation event in time for fitness estimation
- ``nRep`` (_optional_) - number of repeated mutation(if any) starting from 1
- ``theta`` (_optional_) - value of theta (mean cancer cells mutation rate). Default is 0.01.

## Hot to run

In console type and change to actual parameters:

``matlab -nodisplay -nodesktop -r "n=<number>;m=<number>;gv_file='<input_gv_tree>';output='<output_file>';method='<method_name>';nRep=<number>;theta=<number>;SCIFIL"``

Example:

``matlab -nodisplay -nodesktop -r "gv_file='data/dataHou18_map0.gv';names_file='data/dataHou18names.txt';nRep=3;n=58;m=18;SCIFIL"``

Execution result:
![Example](img/example.PNG)

There will be two figures. First represent mutation tree, second fitness landscape.


### Output file

Output contains only one line - calculated fitness of mutations in the same order as in names file or as their numbers in gv file (if names file is not specified). Repeated mutation's fitness will be the last one

```
1.0051 1.0119 1.0206 1.0289 1.0385 1.0490 1.0605 1.0745 1.0904 1.1068 1.1354 1.1107 1.1475 1.1905 1.2466 1.3785 
```