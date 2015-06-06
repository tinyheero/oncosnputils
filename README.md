# oncosnputils
An R package for post-processing the results of OncoSNP

To install this R package, you will need [devtools](http://cran.r-project.org/web/packages/devtools/index.html). Then open R and type:

```r
install_github("tinyheero/oncosnputils")
```

# Key Functionality

The oncosnputils R package main purposes is to supplement some of the output data from OncoSNP. 

For instance, the function `AddLRRBAF2OncosnpCNV` functions will add Log R Ratio (LRR) and B-allele frequeny (BAF) information to each CNV segment in the OncoSNP \*.cnvs file. Additionally, the `AddOncosnp2PennCNVProbe` adds OncoSNP state information to the raw probe data from the PennCNV-Affy protocol. This is useful for plotting the LRR and BAF plots when you want to add colors to the probes.

Also, several loading functions are provided:

* LoadOncoSnpQcFile
* LoadOncosnpCnvFile
* LoadPennCNVProbeData

These provide quick and easy loading (thanks to readr) of these tabular OncoSNP output files and PennCNV probe data. Also, the column names will be renamed to help faciliate downstream R analyses.
