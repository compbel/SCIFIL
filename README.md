# SCIFIL
Single Cell Inference of FItness Landscape

We propose a computational method for in vivo inference of clonal
selection and estimate of fitness landscapes of heterogeneous cancer cell populations from
single cell sequencing data.

It takes single cell data, mutation tree and estimates finesses of all mutations.

## Parameters
- ``matrix_file`` - path to file with cells mutation profile
- ``gv_file``  - path to file with mutation tree
- ``output`` - path for output file
- ``method`` - method to use for fitness calculation:
    - ``heuristic`` to use heuristic
    - ``brute_force`` to use brute force in finding order of mutation event in time for fitness estimation

## Hot to run

In console type:

``matlab -nodisplay -nodesktop -r "matrix_file='<input_mutations_csv>';gv_file='<input_gv_tree>';output='test.out';method='heuristic';fit_from_edges_list"``

Example:

``matlab -nodisplay -nodesktop -r "matrix_file='data/scite/test1_30.csv';gv_file='data/scite/test1_30.gv';output='test.out';method='heuristic';fit_from_edges_list"``

### Input files
There are two input files:
- csv with mutations where each row represents mutation, each column represents cell
- gv file with mutation tree
Examples can be found in data/scite/ folder


### Output file

Output contains only one line - calculated fitness of mutations in the same order as in input csv file:

```
1.0051 1.0119 1.0206 1.0289 1.0385 1.0490 1.0605 1.0745 1.0904 1.1068 1.1354 1.1107 1.1475 1.1905 1.2466 1.3785 
```