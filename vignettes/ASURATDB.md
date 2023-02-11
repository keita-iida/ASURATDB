ASURATDB
================
Keita Iida
2023-02-11

-   <a href="#1-installations" id="toc-1-installations">1 Installations</a>
-   <a href="#2-collect-disease-ontology-database-for-human-cells"
    id="toc-2-collect-disease-ontology-database-for-human-cells">2 Collect
    Disease Ontology database for human cells</a>
-   <a href="#3-collect-cell-ontology-database-for-human-and-mouse-cells"
    id="toc-3-collect-cell-ontology-database-for-human-and-mouse-cells">3
    Collect Cell Ontology database for human and mouse cells</a>
-   <a href="#4-collect-gene-ontology-database-for-human-and-mouse-cells"
    id="toc-4-collect-gene-ontology-database-for-human-and-mouse-cells">4
    Collect Gene Ontology database for human and mouse cells</a>
-   <a href="#5-collect-kegg-database-for-human-and-mouse-cells"
    id="toc-5-collect-kegg-database-for-human-and-mouse-cells">5 Collect
    KEGG database for human and mouse cells</a>
-   <a href="#6-collect-msigdb-for-human-cells"
    id="toc-6-collect-msigdb-for-human-cells">6 Collect MSigDB for human
    cells</a>
    -   <a href="#61-h-hallmark-gene-sets" id="toc-61-h-hallmark-gene-sets">6.1
        H: hallmark gene sets</a>
    -   <a href="#62-c3-regulatory-target-gene-sets-tftgtrd"
        id="toc-62-c3-regulatory-target-gene-sets-tftgtrd">6.2 C3: regulatory
        target gene sets (TFT:GTRD)</a>
    -   <a href="#63-cell-types-in-msigdb" id="toc-63-cell-types-in-msigdb">6.3
        Cell types in MSigDB</a>
-   <a href="#7-collect-cellmarker-for-human-cells"
    id="toc-7-collect-cellmarker-for-human-cells">7 Collect CellMarker for
    human cells</a>
-   <a href="#8-create-a-custom-built-database"
    id="toc-8-create-a-custom-built-database">8 Create a custom-built
    database</a>
    -   <a href="#81-combine-co-and-msigdb-for-human-cells"
        id="toc-81-combine-co-and-msigdb-for-human-cells">8.1 Combine CO and
        MSigDB for human cells</a>
    -   <a href="#82-combine-co-msigdb-and-cellmarker-for-human-cells"
        id="toc-82-combine-co-msigdb-and-cellmarker-for-human-cells">8.2 Combine
        CO, MSigDB, and CellMarker for human cells</a>
    -   <a href="#83-combine-do-co-and-msigdb-for-human-cells"
        id="toc-83-combine-do-co-and-msigdb-for-human-cells">8.3 Combine DO, CO,
        and MSigDB for human cells</a>
    -   <a href="#84-combine-do-co-msigdb-and-cellmarker-for-human-cells"
        id="toc-84-combine-do-co-msigdb-and-cellmarker-for-human-cells">8.4
        Combine DO, CO, MSigDB, and CellMarker for human cells</a>
-   <a href="#9-session-information" id="toc-9-session-information">9
    Session information</a>

# 1 Installations

Attach necessary libraries:

``` r
library(ASURATDB)
```

<br>

# 2 Collect Disease Ontology database for human cells

``` r
library(DOSE) # For using `data(DO2EG)`
```

ASURATDB function `format_DO()` reformats a Disease Ontology database.

``` r
data(DO2EG)
dict_DO <- enrichDO(unlist(DO2EG), ont = "DO", pvalueCutoff = 1,
                    pAdjustMethod = "BH", minGSSize = 0, maxGSSize = 1e+10,
                    qvalueCutoff = 1, readable = FALSE)
human_DO <- format_DO(dict = dict_DO@result, all_geneIDs = dict_DO@gene,
                      orgdb = org.Hs.eg.db::org.Hs.eg.db)
# Save data.
# save(human_DO, file = "genes2bioterm/20201213_human_DO.rda")
```

The data were stored in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 3 Collect Cell Ontology database for human and mouse cells

ASURATDB functions `collect_CO()` and `format_CO()` load a Cell Ontology
database using ontoProc package and reformat the database, respectively.

**Tips:** As of December 2020, Cell Ontology database might not be
complete enough for some biological contexts. For example, well-known
marker genes for pancreatic beta cell, Ins1 and Ins2, were not
registered for “type B pancreatic cell” with ID “CL:0000169”.

``` r
# Human
dict_CO <- collect_CO(orgdb = org.Hs.eg.db::org.Hs.eg.db)
human_CO <- format_CO(dict = dict_CO, orgdb = org.Hs.eg.db::org.Hs.eg.db)
# Save data.
# save(human_CO, file = "genes2bioterm/20201213_human_CO.rda")

# Mouse
dict_CO <- collect_CO(orgdb = org.Mm.eg.db::org.Mm.eg.db)
mouse_CO <- format_CO(dict = dict_CO, orgdb = org.Mm.eg.db::org.Mm.eg.db)
# Save data.
# save(mouse_CO, file = "genes2bioterm/20201211_mouse_CO.rda")
```

The data were stored in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 4 Collect Gene Ontology database for human and mouse cells

ASURATDB functions `collect_GO()` and `format_GO()` load a Gene Ontology
database using clusterProfiler package and reformat the database,
respectively. Currently, only human and mouse data are acceptable.

``` r
# Human
dict_GO <- collect_GO(orgdb = org.Hs.eg.db::org.Hs.eg.db)
human_GO <- format_GO(dict = dict_GO, orgdb = org.Hs.eg.db::org.Hs.eg.db)
# Human reduced
human_GO_red <- human_GO
onts <- c("MF", "BP", "CC")
for(i in seq_along(onts)){
  ids <- human_GO[[onts[i]]][which(human_GO[[onts[i]]]$Count >= 2), ]$ID
  mat <- human_GO$similarity_matrix[[onts[i]]][ids, ids]
  human_GO_red$similarity_matrix[[onts[i]]] <- mat
}
# Save data.
# save(human_GO_red, file = "genes2bioterm/20201213_human_GO_red.rda")

# Mouse
dict_GO <- collect_GO(orgdb = org.Mm.eg.db::org.Mm.eg.db)
mouse_GO <- format_GO(dict = dict_GO, orgdb = org.Mm.eg.db::org.Mm.eg.db)
# Mouse reduced
mouse_GO_red <- mouse_GO
onts <- c("MF", "BP", "CC")
for(i in seq_along(onts)){
  ids <- mouse_GO[[onts[i]]][which(mouse_GO[[onts[i]]]$Count >= 2), ]$ID
  mat <- mouse_GO$similarity_matrix[[onts[i]]][ids, ids]
  mouse_GO_red$similarity_matrix[[onts[i]]] <- mat
}
# Save data.
# save(mouse_GO_red, file = "genes2bioterm/20201211_mouse_GO_red.rda")
```

The data were stored in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 5 Collect KEGG database for human and mouse cells

ASURATDB functions `collect_KEGG()` and `format_KEGG()` load a KEGG
database using KEGGREST package via the internet and reformat the
database, respectively.

The arguments of `collect_KEGG()` are `organism` and `categories`. Here,
`organism` must obey the naming rule of
[KEGG](http://rest.kegg.jp/list/organism) (see `KEGGREST` function
`listDatabases()`) and `categories` must be one of `"pathway"`,
`"module"`, and `"drug"` (only for human) in the current version.

``` r
# Human
dict_KEGG <- collect_KEGG(organism = "hsa", categories = c("pathway"))
human_KEGG <- format_KEGG(dict = list(pathway = dict_KEGG[["pathway"]][["success"]]),
                          orgdb = org.Hs.eg.db::org.Hs.eg.db)
# Save data.
# save(human_KEGG, file = "genes2bioterm/20201213_human_KEGG.rda")

# Mouse
dict_KEGG <- collect_KEGG(organism = "mmu", categories = c("pathway"))
mouse_KEGG <- format_KEGG(dict = list(pathway = dict_KEGG[["pathway"]][["success"]]),
                          orgdb = org.Mm.eg.db::org.Mm.eg.db)
# Save data.
# save(mouse_KEGG, file = "genes2bioterm/20201211_mouse_KEGG.rda")

# Human (drug)
dict_KEGG_drug <- collect_KEGG(organism = "hsa", categories = c("drug"))
human_KEGG_drug <- format_KEGG(dict = list(drug = dict_KEGG_drug[["drug"]][["success"]]),
                               orgdb = org.Hs.eg.db::org.Hs.eg.db)
# Save data.
# save(human_KEGG_drug, file = "genes2bioterm/20221102_human_KEGG_drug.rda")
```

Note `collect_KEGG()` uses `KEGGREST` function `keggGet()`, which may
produce both successful and unsuccessful results. The data were stored
in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 6 Collect MSigDB for human cells

## 6.1 H: hallmark gene sets

Load databases, where category is “H” (hallmark gene sets) and species
is human (cf. `msigdbr::msigdbr_species()`).

``` r
dbtable <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
```

Reformat the database.

``` r
dbtable_gsetID <- dbtable[, which(colnames(dbtable) %in% c("gs_name", "gs_id"))]
dbtable_gsetID <- unique(dbtable_gsetID)
dbtable_geneID <- split(x = dbtable$human_entrez_gene, f = dbtable$gs_name)
dbtable_symbol <- split(x = dbtable$gene_symbol, f = dbtable$gs_name)
stopifnot(identical(length(dbtable_geneID), length(dbtable_symbol)))

res <- c("ID", "Description", "Count", "Gene", "GeneID", "IC")
res <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, res)))
for(i in 1:length(dbtable_geneID)){
  res <- rbind(res, data.frame(
    ID = dbtable_gsetID$gs_id[i],
    Description = dbtable_gsetID$gs_name[i],
    IC = NA,
    Count = length(dbtable_geneID[[i]]),
    Gene = paste(dbtable_symbol[[i]], collapse = "/"),
    GeneID = paste(dbtable_geneID[[i]], collapse = "/")))
}
human_MSigDB_Hallmark <- list(hallmark = res)
# Save data.
# save(human_MSigDB_Hallmark, file = "genes2bioterm/20230127_human_MSigDB_Hallmark.rda")
```

The data were stored in the following repositories:

-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

## 6.2 C3: regulatory target gene sets (TFT:GTRD)

Load databases, where category is “C3” (regulatory target gene sets) and
species is human (cf. `msigdbr::msigdbr_species()`).

``` r
dbtable <- msigdbr::msigdbr(species = "Homo sapiens", category = "C3")
dbtable <- dbtable[which(dbtable$gs_subcat == "TFT:GTRD"), ]
```

Reformat the database.

``` r
dbtable_gsetID <- dbtable[, which(colnames(dbtable) %in% c("gs_name", "gs_id"))]
dbtable_gsetID <- unique(dbtable_gsetID)
dbtable_geneID <- split(x = dbtable$human_entrez_gene, f = dbtable$gs_name)
dbtable_symbol <- split(x = dbtable$gene_symbol, f = dbtable$gs_name)
stopifnot(identical(length(dbtable_geneID), length(dbtable_symbol)))

res <- c("ID", "Description", "Count", "Gene", "GeneID", "IC")
res <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, res)))
for(i in 1:length(dbtable_geneID)){
  res <- rbind(res, data.frame(
    ID = dbtable_gsetID$gs_id[i],
    Description = dbtable_gsetID$gs_name[i],
    IC = NA,
    Count = length(dbtable_geneID[[i]]),
    Gene = paste(dbtable_symbol[[i]], collapse = "/"),
    GeneID = paste(dbtable_geneID[[i]], collapse = "/")))
}
human_MSigDB_GTRD <- list(GTRD = res)
# Save data.
# save(human_MSigDB_GTRD, file = "genes2bioterm/20230211_human_MSigDB_GTRD.rda")
```

The data were stored in the following repositories:

-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

## 6.3 Cell types in MSigDB

Load databases.

``` r
dbtable <- clustermole::clustermole_markers()
```

``` r
sort(unique(dbtable$db))
```

    [1] "ARCHS4"     "CellMarker" "MSigDB"     "PanglaoDB"  "SaVanT"     "TISSUES"   
    [7] "xCell"

Select species and databases.

``` r
dbtable <- dbtable[which(dbtable$species == "Human"), ]
dbtable <- dbtable[which(dbtable$db == "MSigDB"),]
dbtable$geneID <- NA
```

Change gene symbols into entrez IDs.

``` r
dictionary <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                    key = dbtable$gene_original,
                                    columns = c("SYMBOL", "ENTREZID"),
                                    keytype = "SYMBOL")
dictionary <- dictionary[!duplicated(dictionary$SYMBOL), ]
dictionary <- dictionary[which(!is.na(dictionary$SYMBOL)),]
for(i in 1:nrow(dbtable)){
  gene <- dbtable$gene_original[i]
  inds <- which(dictionary$SYMBOL == gene)
  dbtable$geneID[i] <- dictionary[inds,]$ENTREZID
}
```

Reformat the database. Here, the identifier of each biological term are
named “MSigDBID.”

``` r
dbtable_geneID <- split(x = dbtable$geneID, f = dbtable$celltype)
dbtable_symbol <- split(x = dbtable$gene_original, f = dbtable$celltype)
stopifnot(identical(length(dbtable_geneID), length(dbtable_symbol)))

res <- c("ID", "Description", "Count", "Gene", "GeneID", "IC")
res <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, res)))
for(i in 1:length(dbtable_geneID)){
  res <- rbind(res, data.frame(
    ID = paste("MSigDBID:", i, sep = ""),
    Description = names(dbtable_geneID)[i],
    IC = NA,
    Count = length(dbtable_geneID[[i]]),
    Gene = paste(dbtable_symbol[[i]], collapse = "/"),
    GeneID = paste(dbtable_geneID[[i]], collapse = "/")))
}
human_MSigDB <- list(cell = res)
# Save data.
# save(human_MSigDB, file = "genes2bioterm/20220308_human_MSigDB.rda")
```

The data were stored in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 7 Collect CellMarker for human cells

Load databases.

``` r
dbtable <- clustermole::clustermole_markers()
```

``` r
sort(unique(dbtable$db))
```

    [1] "ARCHS4"     "CellMarker" "MSigDB"     "PanglaoDB"  "SaVanT"     "TISSUES"   
    [7] "xCell"

Select species and databases.

``` r
dbtable <- dbtable[which(dbtable$species == "Human"), ]
dbtable <- dbtable[which(dbtable$db == "CellMarker"),]
dbtable$geneID <- NA
```

Change gene symbols into entrez IDs.

``` r
dictionary <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                    key = dbtable$gene_original,
                                    columns = c("SYMBOL", "ENTREZID"),
                                    keytype = "SYMBOL")
dictionary <- dictionary[!duplicated(dictionary$SYMBOL), ]
dictionary <- dictionary[which(!is.na(dictionary$SYMBOL)),]
for(i in 1:nrow(dbtable)){
  gene <- dbtable$gene_original[i]
  inds <- which(dictionary$SYMBOL == gene)
  dbtable$geneID[i] <- dictionary[inds,]$ENTREZID
}
```

Reformat the database. Here, the identifier of each biological term are
named “CellMarkerID.”

``` r
dbtable_geneID <- split(x = dbtable$geneID, f = dbtable$celltype)
dbtable_symbol <- split(x = dbtable$gene_original, f = dbtable$celltype)
stopifnot(identical(length(dbtable_geneID), length(dbtable_symbol)))

res <- c("ID", "Description", "Count", "Gene", "GeneID", "IC")
res <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, res)))
for(i in 1:length(dbtable_geneID)){
  res <- rbind(res, data.frame(
    ID = paste("CellMarkerID:", i, sep = ""),
    Description = names(dbtable_geneID)[i],
    IC = NA,
    Count = length(dbtable_geneID[[i]]),
    Gene = paste(dbtable_symbol[[i]], collapse = "/"),
    GeneID = paste(dbtable_geneID[[i]], collapse = "/")))
}
human_CellMarker <- list(cell = res)
# Save data.
# save(human_CellMarker, file = "genes2bioterm/20220308_human_CellMarker.rda")
```

The data were stored in the following repositories:

-   [DOI:10.6084/m9.figshare.19102598](https://figshare.com/s/0599d2de970c2deb675c)
-   [Github ASURATDB](https://github.com/keita-iida/ASURATDB)

<br>

# 8 Create a custom-built database

## 8.1 Combine CO and MSigDB for human cells

Create a cell type-related database by combining Cell ontology and
MSigDB databases for analyzing human single-cell transcriptome data.

``` r
urlpath <- "https://github.com/keita-iida/ASURATDB/blob/main/genes2bioterm/"
load(url(paste0(urlpath, "20201213_human_CO.rda?raw=true")))
load(url(paste0(urlpath, "20220308_human_MSigDB.rda?raw=true")))
res <- rbind(human_CO[["cell"]], human_MSigDB[["cell"]])
human_CB <- list(cell = res)
```

<br>

## 8.2 Combine CO, MSigDB, and CellMarker for human cells

Create a cell type-related database by combining Cell ontology, MSigDB,
and CellMarker databases for analyzing human single-cell transcriptome
data.

``` r
urlpath <- "https://github.com/keita-iida/ASURATDB/blob/main/genes2bioterm/"
load(url(paste0(urlpath, "20201213_human_CO.rda?raw=true")))
load(url(paste0(urlpath, "20220308_human_MSigDB.rda?raw=true")))
load(url(paste0(urlpath, "20220304_human_CellMarker.rda?raw=true")))
res <- do.call("rbind", list(human_CO[["cell"]], human_MSigDB[["cell"]],
                             human_CellMarker[["cell"]]))
human_CB <- list(cell = res)
```

<br>

## 8.3 Combine DO, CO, and MSigDB for human cells

Create a cell type-related database by combining Disease Ontology, Cell
ontology and MSigDB databases for analyzing complex human single-cell
transcriptome data.

``` r
urlpath <- "https://github.com/keita-iida/ASURATDB/blob/main/genes2bioterm/"
load(url(paste0(urlpath, "20201213_human_DO.rda?raw=true")))
load(url(paste0(urlpath, "20201213_human_CO.rda?raw=true")))
load(url(paste0(urlpath, "20220308_human_MSigDB.rda?raw=true")))
res <- do.call("rbind", list(human_DO[["disease"]], human_CO[["cell"]],
                             human_MSigDB[["cell"]]))
human_CB <- list(cell = res)
```

<br>

## 8.4 Combine DO, CO, MSigDB, and CellMarker for human cells

Create a cell type-related database by combining Disease Ontology, Cell
ontology, MSigDB, and CellMarker databases for analyzing complex human
single-cell transcriptome data.

``` r
urlpath <- "https://github.com/keita-iida/ASURATDB/blob/main/genes2bioterm/"
load(url(paste0(urlpath, "20201213_human_DO.rda?raw=true")))
load(url(paste0(urlpath, "20201213_human_CO.rda?raw=true")))
load(url(paste0(urlpath, "20220308_human_MSigDB.rda?raw=true")))
load(url(paste0(urlpath, "20220304_human_CellMarker.rda?raw=true")))
res <- do.call("rbind", list(human_DO[["disease"]], human_CO[["cell"]],
                             human_MSigDB[["cell"]], human_CellMarker[["cell"]]))
human_CB <- list(cell = res)
```

<br>

# 9 Session information

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur ... 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.2.1  fastmap_1.1.0   cli_3.6.0       tools_4.2.1    
#>  [5] htmltools_0.5.4 rstudioapi_0.14 yaml_2.3.7      rmarkdown_2.20 
#>  [9] knitr_1.42      xfun_0.37       digest_0.6.31   rlang_1.0.6    
#> [13] evaluate_0.20
```
