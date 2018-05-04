# bambi
R package for querying 10x WGS and single-cell BAMs

## Dependencies


```{r}
devtools::install_github('mskilab/gUtils')
```

```{r}
devtools::install_github('mskilab/bamUtils')
```


cbambi
https://github.com/mskilab/cbambi

## bambi commands


Instantiate a bambi object:

```{r}
hcc1143_subset = bambi$new(bam_file = "subsetHCC1143_phased_possorted0001.bam", bamdb_path="subsetHCC1143_phased_possorted0001_lmdb")
```

Methods:

* grab_bx()
    * grab_bx(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)

* grab_cb()
    * grab_cb(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)

* grab_ub()
    * grab_ub(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)

* fetch_by_tag()
    * fetch_by_tag(tag, tag_queries, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)

* fetch_by_multiple_tag()? 
    * remove
    * the only downside is that `fetch_by_tag()` might take a performance hit....





