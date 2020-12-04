This folder contains the data and code needed to reproduce the language-change example (Sec. 2.4). 

- Fig5-ling-example.r: the code for running Clarity on the language-change data
- folder "samples": the similarity matrices serving as the input to Clarity; each file contains two matrices, one for phonetic, the other for lexical similarity. 
- header-NorthEuraLex.tsv: header line for these data matrices

The similarity matrices were computed using the methods described by Dellert (2018), with the code published here under GPL-3.0: https://github.com/jdellert/iwsa  Raw data are from the open resource NorthEuraLex v.0.90, available at http://northeuralex.org/  The methods are explained in some detail in Sec. 4.7 of the paper; consult the original publication Dellert (2018) for more information. 