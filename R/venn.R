#' Make two-circle Venn diagram
#'
#' @param set1 a character vector of values in the 1st set. Required
#' @param set2 a character vector of values in the 2nd set. Required
#' @param names a character vector of length two containing names of the sets. Default: c("set 1", "set 2")
#' @param title title of the Venn diagram. Default: "Venn diagram".
#' 
#' @return a data frame with three columns used to get overlap counts
#' @export
#' @examples
#' \dontrun{
#' #define some sets
#' e <- c(1, 2, 3, 4, 5)
#' f <- c(3, 4, 5, 6, 7)
#' Venn2(e, f, c("e", "f"), "E and F")
#' }
# 

Venn2<-function(set1, set2, names = c("set 1", "set 2"), title = "Venn diagram")
{
	stopifnot(length(names) == 2)
	#Form universe as union of sets
	universe<-sort(unique(c(set1,set2)))
	Counts<-matrix(0,nrow=length(universe),ncol=2)
	colnames(Counts)<-names
	for(i in 1:length(universe))
	{
		Counts[i,1]<- universe[i] %in% set1
		Counts[i,2]<- universe[i] %in% set2
	}
	par(mar=c(1,1,2,1),oma=c(1,1,2,1))
	limma::vennDiagram(limma::vennCounts(Counts));
	mtext(title,outer=T,line=1);
	tab<-data.frame(universe,Counts,stringsAsFactors=FALSE)
	colnames(tab)<-c("id",names);
	return(tab);
}

#' Make three-circle Venn diagram
#'
#' @param set1 a character vector of values in the 1st set. Required
#' @param set2 a character vector of values in the 2nd set. Required
#' @param set3 a character vector of values in the 3rd set. Required
#' @param names a character vector of length two containing names of the sets. Default: c("set 1", "set 2", "set 3")
#' @param title title of the Venn diagram. Default: "Venn diagram".
#' 
#' @return a data frame with four columns used to get overlap counts
#' @export
#' @examples
#' \dontrun{
#' #define some sets
#' e <- c(1, 2, 3, 4, 5)
#' f <- c(3, 4, 5, 6, 7)
#' g <- c(1, 7, 8)
#' Venn3(e, f, g, c("e", "f", "g"), "E, F and G")
#' }
#
Venn3 <- function(set1, set2, set3, names = c("set 1", "set 2", "set 3"), title = "Venn diagram")
{
	stopifnot( length(names) == 3)
	# Form universe as union of all three sets
	universe <- sort( unique( c(set1, set2, set3) ) )
	Counts <- matrix(0, nrow=length(universe), ncol=3)
	colnames(Counts) <- names

	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
		Counts[i,3] <- universe[i] %in% set3
	}

	par(mar=c(1,1,2,1),oma=c(1,1,2,1))
	limma::vennDiagram( limma::vennCounts(Counts) )
	mtext(title,outer=T,line=1);
	tab<-data.frame(universe,Counts,stringsAsFactors=FALSE)
	colnames(tab)<-c("id",names);
	return(tab);
}


