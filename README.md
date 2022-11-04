# MARS
Under the high-dimensional setting, we proposed a second-order algorithm (MARS) for estimating sparse precision matrices. Note that our algorithm is designed for the high dimensional setting, where the dimension $p$ is much larger than the sample size $n$.

# References
Qian LI, Binyan Jiang, and Defeng Sun. ["MARS: A second-order reduction algorithm for high-dimensional sparse precision matrices estimation"](https://www.polyu.edu.hk/ama/profile/dfsun/MARS0624.pdf). 

# Getting started
These instructions will give you a toy example for implementing the package.

## Preparation
Before proceeding to use the package, you need install Rtools under Windows and Xcode under MacOS for compiling the c++ codes and also some necessary R packages.

```
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("Matrix")
```

## Installing MARS

```
library("devtools")
devtools::install_github("QianLI-QL/MARS")
```

## Example

After installing MARS, you can test the insulation by typing the following into R.
```
library('Rcpp')
library('RcppArmadillo')
library('Matrix')
library('MARS')
solution <- MARS(matrix(runif(3000*50), 3000, 50), stopmethod = "bigs", Lambdapath = seq(1, 0.5, -0.01), printmain = TRUE)
```
