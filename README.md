# SCIFIL
Single Cell Inference of FItness Landscape

We propose a computational method for in vivo inference of clonal
selection and estimate of fitness landscapes of heterogeneous cancer cell populations from
single cell sequencing data.

## Hot to run

In console type:

``matlab -nodisplay -nodesktop -r "input='<input_file>';output='<output_file>'';fit_from_edges_list"``

Example:

``matlab -nodisplay -nodesktop -r "input='input.test';output='test.out';fit_from_edges_list"``

### Input format
Input file has following structure:
```
<number of mutations> <wild_type_frequency>
<mutation_id_1> <parent_1> <frequency_1>
...
<mutation_id_n> <parent_n> <frequency_n>
```

Example:
```
16 112
1 0 4
2 1 1
3 2 2
4 3 1
5 4 1
6 5 1
7 6 1
8 7 3
9 8 10
10 9 10
11 10 9
12 10 0.001
13 12 1
14 13 4
15 14 3
16 15 16
```
16 112 - 16 mutations, frequency of wild type is 112
1 0 4 - mutation 1 is a child of wild type(1) with frequency 4
2 1 1 - mutation 2 is a child of mutation 1 with frequency 1
...

This example corresponds to following tree:

![Example](img/example.PNG)

### Output format

There are two lines:
- order of mutations occurrences in time
- fitness of mutations in initial order (i.e. first number - fitness of first mutation in input, 
n-th number - fitness of n-th mutation in input)

```
1 2 3 4 5 6 7 8 9 10 12 13 11 14 15 16 
1.0051 1.0119 1.0206 1.0289 1.0385 1.0490 1.0605 1.0745 1.0904 1.1068 1.1354 1.1107 1.1475 1.1905 1.2466 1.3785 
```