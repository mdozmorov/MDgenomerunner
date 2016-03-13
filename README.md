MD helper functions
===

Package created with the help of [Dave Tang's tutorial](http://davetang.org/muse/2015/02/04/bed-granges/#more-5066).

To install:

```r
library(devtools)
install_github('mdozmorov/MDmisc')
library(MDmisc)
```

# Misc functions

`gene_enrichment` - gene set enrichment analysis on a list of gene names or EntrezIDs. Currently, Homo Sapiens only.

`unembed` - extracts multiple values embedded in rows. E.g. one row with ABC11 /// BCD22 variable will be split into two separate entries, creating two ABC11 and BCD22 rows with other values equal to the original row.

`get_pubmed_graph` - searches a term or phase in PubMed within year limits, and plot a barplot of counts.

`TCGA_corr` - retrieves correlation statistics between RSEM gene expression of a gene of interest and all other genes in a cancer sybtype from TCGA.

`Venn2`, `Venn3`, `Venn5` - plot 2, 3, 5-circle Venn diagrams.

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


# Misc notes

### Rmd header

	```{r setup, echo=FALSE, message=FALSE, warning=FALSE}
	# Set up the environment
	library(knitr)
	opts_chunk$set(cache.path='cache/', fig.path='img/', cache=F, tidy=T, fig.keep='high', echo=F, dpi=100, warnings=F, message=F, comment=NA, warning=F, results='as.is') #out.width=700, 
	library(pander)
	panderOptions('table.split.table', Inf)
	set.seed(1)
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