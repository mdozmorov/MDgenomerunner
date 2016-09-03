MD helper functions
===

Package created with the help of [Dave Tang's tutorial](http://davetang.org/muse/2015/02/04/bed-granges/#more-5066).

To install:

```r
library(devtools)
install_github('mdozmorov/MDmisc')
library(MDmisc)
```

Or, from command line:

```bash
sudo R -e "devtools::install_github('mdozmorov/MDmisc')"
```

# Misc functions

`gene_filter` - filters expression matrix by two genefilter criteria

`gene_enrichment` - gene set enrichment analysis on a list of gene names or EntrezIDs. Currently, Homo Sapiens and Mus Musculus only. Now, supporting GO, KEGG, and [msigdf](https://github.com/stephenturner/msigdf) enrichments.

`unembed` - extracts multiple values embedded in rows. E.g. one row with ABC11 /// BCD22 variable will be split into two separate entries, creating two ABC11 and BCD22 rows with other values equal to the original row.

`get_pubmed_graph` - searches a term or phase in PubMed within year limits, and plot a barplot of counts.

`TCGA_corr` - retrieves correlation statistics between RSEM gene expression of a gene of interest and all other genes in a cancer sybtype from TCGA.

`Venn2`, `Venn3`, `Venn5` - plot 2, 3, 5-circle Venn diagrams.

`mtx_rand` - randomizes a numerical matrix using several techniques

`t.test2` - t-test from a summary statistics

### GenomeRunner-specific main functions

`gr_load_data` - loads enrichment analysis matrix(es) and remove non-informative enrichment results.

`gr_clusters` - defines clusters by cutting a dendrogram by specified height.

`gr_degfs` - compares the enrichment values between pairs of clusters using Wilcoxon (aka Mann-Whitney) test.

`gf_cellspecific` - performs cell type-specific enrichment analysis.

### GenomeRunner-specific ancillary functions

`gr_transform` -  converts a matrix of p-values (with "-" indicating depletion) into a -log10-transformed matrix with sign preserved.

`gr_untransform` -  converts a matrix of -log10-transformed p-values (with "-" indicating depletion) into a matrix of regular p-values with sign preserved.

`gr_trimnames` - trims long column/row names to a defined length.

`gr_maxmin` - extracts max/min correlation of each object (column) with other objects. Each column in an n x m matrix is correlated with other columns, and the names of the columns maximally/minimally correlated with the target column are assembled.

`gr_promoter_extract` - given a vector of gene EntrezIDs, or symbols, extract genomic coordinates of their promoters. Promoters can be defined as (by default) regions 2,000 bp upstream and 500 bp downstream of gene' transcription start site. Currently, hg19 only.

`gr_gene_extract` - given a vector of gene EntrezIDs, or symbols, extract genomic coordinates of the full genes. Currently, hg19 only.

`gr_plot` - master function to visualize enrichment results outside the server environment.

`gr_getV1OddsRatioPvalMatrix` - converts a log file ofthe results obtained with GenomeRunner V1 (standalone) to "matrix_PVAL/OR.txt" files compatible with GenomeRunner V2 (web)

# Misc notes

### Global settings

`options(stringsAsFactors = FALSE)` - place at the beginning of R code

### YAML header

[R Markdown v2 formatting options](http://rmarkdown.rstudio.com/html_document_format.html#overview)

	---
	title: "Demo Document"
	output:
	  html_document:
	    toc: true
	    theme: united
	date: "`r Sys.Date()`"
	author: "Author's Name"
	---

### Rmd header

	```{r setup, echo=FALSE, message=FALSE, warning=FALSE}
	# Set up the environment
	library(knitr)
	opts_chunk$set(cache.path='cache/', fig.path='img/', cache=F, tidy=T, fig.keep='high', echo=F, dpi=100, warnings=F, message=F, comment=NA, warning=F, results='as.is') #out.width=700, 
	library(pander)
	panderOptions('table.split.table', Inf)
	set.seed(1)
	library(dplyr)
	options(stringsAsFactors = FALSE)
	```

### Rmd footer

  	```{r session_info}
  	diagnostics <- devtools::session_info()
    platform <- data.frame(diagnostics$platform %>% unlist, stringsAsFactors = FALSE)
    colnames(platform) <- c("description")
    pander(platform)
    
    packages <- as.data.frame(diagnostics$packages)
    pander(packages[ packages$`*` == "*", ])
  	```

### Tidy up code

`formatR::tidy_app()`


# Upgrading R & Bioconductor on Mac OS X

When upgrading to the latest X.Y.Z R vestion, it is possible to rename the `/Library/Frameworks/R.framework/Versions/X.Y` folder. This will make all packages available to the latest R version, but often causes unforeseen incompatibilities, problems with upgrading Bioconductor, and other unpredictable behaviour. A more laborous manual upgrade may be simpler and cleaner.

1. Rename the current R library folder, e.g., `/Library/Frameworks/R.framework/Versions/3.2` to `/Library/Frameworks/R.framework/Versions/3.2.old`
2. Install the latest version of [R](https://www.r-project.org/)
3. Install the latest version of [Bioconductor](https://www.bioconductor.org/install/)
4. Go to the terminal, and note the new R folder, `d1=/Library/Frameworks/R.framework/Versions/3.3/Resources/library`
5. Note the old R folder, `d2=/Library/Frameworks/R.framework/Versions/3.2.old/Resources/library`
6. Get the list of R packages you need to install from the previous installation, `comm -13 <(find $d1 -type d -maxdepth 1 -exec basename {} \; | sort) <(find $d2 -type d -maxdepth 1 -exec basename {} \; | sort) > ~/Desktop/R_packages_to_upgrade.txt`
7. Use the `R_packages_to_upgrade.txt` list to install necessary R packages into the fresh R/Bioconductor installation. The majority of them are dependencies, so use intelligence to install what you need.

	```{r}
	source("https://bioconductor.org/biocLite.R")

	library(BiocInstaller)

	biocLite(c("ChIPseeker", "clusterProfiler", "org.Hs.eg.db", "Rgraphviz", "pathview", "genefilter", "Category", "edgeR", "sva", "ReactomePA", "GOstats", "KEGG.db", "methylumi", "lumi", "wateRmelon", "betareg", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", "hgu133a.db", "hgu133a2.db", "impute", "sscore", "snpStats", "bladderbatch", "simpleaffy", "fpc", "GEOquery", "e1071", "pathview", "biclust", "eisa", "ExpressionView", "samr", "WGCNA", "RISmed"))

	biocLite('CellMix', siteRepos = 'http://web.cbio.uct.ac.za/~renaud/CRAN', type='both')

	# brew install openssl
	install.packages(c("devtools", "TCGA2STAT", "pander", "xlsx", "ggrepel", "shiny", "shinyBS", "devtools", "roxygen2", "caret", "kernlab", "pROC", "openxlsx", "XLConnect", "pheatmap", "scatterplot3d", "tsne", "Rtsne", "DT", "pvclust", "Hmisc", "dynamicTreeCut", "apcluster", "rgl", "calibrate"))

	devtools::install_github("mdozmorov/MDmisc")
	devtools::install_github("mdozmorov/annotables")

	install_github('ramnathv/slidify')
	install_github('ramnathv/slidifyLibraries')

	# Install Java, then
	install.packages(c("rJava", "xlsx"))
	```

