# MAIC
``Minimization of Akaike's Information Criterion
in Linear Regression Analysis via Mixed Intger Nonlinear Program"
by Keiji Kimura and Hayato Waki

## Description
The name of this software is MAIC which is written by C language.
MAIC minimizes Akaike's Information Criterion (AIC) in Linear Regression.
cblas.h, clapack.h and SCIP(http://scip.zib.de/) necessary to use MAIC.
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

## Licence
coming soon ....
