# oncosnputils - An R package for post-processing the results of OncoSNP

To install this R package, you will need [devtools](http://cran.r-project.org/web/packages/devtools/index.html). Then open R and type:

```r
install_github("tinyheero/oncosnputils")
```

# Key Functionality

The oncosnputils R package main purposes is to supplement some of the output data from OncoSNP. 

For instance, the function `add_LRR_BAF_to_oncosnp_cnv` functions will add Log R Ratio (LRR) and B-allele frequeny (BAF) information to each CNV segment in the OncoSNP \*.cnvs file. Additionally, the `add_oncosnp_to_penncnv_probe` adds OncoSNP state information to the raw probe data from the PennCNV-Affy protocol. This is useful for plotting the LRR and BAF plots when you want to add colors to the probes.

Also, several loading functions are provided:

* load_oncosnp_qc_file
* load_oncosnp_cnv_file
* load_penncnv_probe_data

These provide quick and easy loading (thanks to data.table::fread) of these tabular OncoSNP output files and PennCNV probe data. Also, the column names will be renamed to help faciliate downstream R analyses.
