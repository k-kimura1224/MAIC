# MAIC

## Description
The name of this software is MAIC, which minimizes Akaike's
Information Criterion (AIC) in Linear Regression. We implement a
relaxator and branching rules by C language. These are called from
SCIP (http://scip.zib.de/), which is a mathematical optimization
software and a branch-and-bound framework. These algorithms are based
on our manuscript ``Minimization of Akaike's Information Criterion in
Linear Regression Analysis via Mixed Intger Nonlinear Program" by
Keiji Kimura and Hayato Waki.

cblas.h, clapack.h and SCIP(http://scip.zib.de/) are required to
install MIAC. See thier web sites for their installation. For
instance, colas.h and clapack.h are available at ubuntu.

All sample datas in this software are obtained from UCI Machine
Learning Repository (http://archive.ics.uci.edu/ml/).


## Compile
Write in Makefile path of the directory of your scip in SCIPDIR
Then, execute ``make" to compile our software.

## Usage
You can solve d1_housing.linereg in ``data" by the following command

 - ./bin/linereg -f data/d1_housing.linereg

## Definition of Data
Our software reads from the second line of data file.
The definition of our data is as follows:
- [2nd line] the number of data
- [3rd line] the number of the explanatory variables
- [4th line] the index of the response variable(Begining is 1)
- [from 5th line] your data

See data files in the directory data for more details.

