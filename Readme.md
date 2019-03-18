# Clarity: An R package for Comparing SimiLARITY Matrices

## Overview

Clarity is an R package for comparing matrices in terms of their *structure*. If the matrices represent dissimilarities between objects, then the matrix can be decomposed into a **structure** and a **relationship** matrix, in the manner of Non-negative Matrix Factorization. Each object is represented in terms of:

* the **Structure**: the topology of a graph (mixture of trees) that generates the objects in a latent space
* a **Relationship** matrix represents the branch lengths between the latent structures.

**Structural comparison** allows for the the structures of two matrices to be compared whilst making minimal assumptions regarding the relationship matrix.

Clarity offers tools to rapidly perform the graph comparison, to identify nodes that are outlying and to identify elements of the matrices that are unusual with respect to one another.

## Installation

Installation should be straightforward with the "devtools" package, which can be installed with `install.packages("devtools")`. Then you can simply call:

```R
library("devtools")
install_github("danjlawson/CLARITY/Clarity")
```

If you would prefer to get the [Clarity tar file](https://github.com/danjlawson/CLARITY/Clarity.tar.gz) you can then install that from a terminal with `R CMD INSTALL Clarity.tar.gz`.

## Usage

You should find usage straightforward if you follow the examples in the Clarity package:
```R
help(Clarity)
```
will point you towards

## Citation

Currently cite the package with `citation("Clarity")`. A paper is under submission.

## License and Authors

`Clarity` and `ClaritySim` are released under the GPLv3 license.

Clarity is written by Daniel Lawson (dan.lawson@bristol.ac.uk).
