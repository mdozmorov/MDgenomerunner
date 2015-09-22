MD helper functions
===

Package created with the help of [Dave Tang's tutorial](http://davetang.org/muse/2015/02/04/bed-granges/#more-5066).

To install:

```r
library(devtools)
install_github('mdozmorov/MDmisc')
library(MDmisc)
```

List of functions:

`gene_enrichment` - gene set enrichment analysis on a list of gene names or EntrezIDs. Currently, Homo Sapiens only.

`unembed` - extracts multiple values embedded in rows. E.g. one row with ABC11 /// BCD22 variable will be split into two separate entries, creating two ABC11 and BCD22 rows with other values equal to the original row.


### GenomeRunner-specific main functions

`gr_load_data` - loads enrichment analysis matrix(es) and remove non-informative enrichment results.

`gr_clusters` - defines clusters by cutting a dendrogram by specified height.

`gr_degfs` - compares the enrichment values between pairs of clusters using Wilcoxon (aka Mann-Whitney) test.

### GenomeRunner-specific ancillary functions

`gr_transform` -  converts a matrix of p-values (with "-" indicating depletion) into a -log10-transformed matrix with sign preserved

`gr_trimnames` - trim long column/row names to a defined length.
