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

#' Make five-circle Venn diagram
#'
#' @param set12345 a character vector of values in the nth set. Required
#' 
#' @return Nothing, just plot the diagram
#' @export
#' @examples
#' \dontrun{
#' # define some sets
#' set1 = letters[1:6];
#' set2 = letters[6:10];
#' set3 = letters[6:15];
#' set4 = letters[6:21];
#' set5 = letters[6:26];
#' Venn5(set1, set2, set3, set4, set5);
#' }
#
Venn5 <- function(set1, set2, set3, set4, set5){
  
  all_ids=unique(c(set1, set2, set3, set4, set5));
  all_n=length(all_ids);
  
  #create counting matrix
  
  all_matrix=cbind(affy=rep(0,all_n),est=rep(0,all_n),glomchip=rep(0,all_n),sage=rep(0,all_n),stanford=rep(0,all_n));
  rownames(all_matrix)=all_ids;
  colnames(all_matrix)=c(substitute(set1), substitute(set2),
                         substitute(set3), substitute(set4), substitute(set5)); 
  all_matrix[all_ids %in% set1, 1]=1;
  all_matrix[all_ids %in% set2, 2]=1;
  all_matrix[all_ids %in% set3, 3]=1;
  all_matrix[all_ids %in% set4, 4]=1;
  all_matrix[all_ids %in% set5, 5]=1;
  
  #elements only in one data set.
  
  only1= apply(all_matrix,1,sum) == 1;
  nA=sum(all_matrix[only1, 1]);
  nB=sum(all_matrix[only1, 2]);
  nC=sum(all_matrix[only1, 3]);
  nD=sum(all_matrix[only1, 4]);
  nE=sum(all_matrix[only1, 5]);
  
  #elements overlapped in two data set.
  
  only2= apply(all_matrix,1,sum) == 2;
  nAB=sum(only2 & all_matrix[,1] & all_matrix[,2]);
  nAC=sum(only2 & all_matrix[,1] & all_matrix[,3]);
  nAD=sum(only2 & all_matrix[,1] & all_matrix[,4]);
  nAE=sum(only2 & all_matrix[,1] & all_matrix[,5]);
  nBC=sum(only2 & all_matrix[,2] & all_matrix[,3]);
  nBD=sum(only2 & all_matrix[,2] & all_matrix[,4]);
  nBE=sum(only2 & all_matrix[,2] & all_matrix[,5]);
  nCD=sum(only2 & all_matrix[,3] & all_matrix[,4]);
  nCE=sum(only2 & all_matrix[,3] & all_matrix[,5]);
  nDE=sum(only2 & all_matrix[,4] & all_matrix[,5]);
  
  #elements overlapped in three data set.
  
  only3= apply(all_matrix,1,sum) == 3;
  nABC=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3]);
  nABD=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,4]);
  nABE=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,5]);
  nACD=sum(only3 & all_matrix[,1] & all_matrix[,3] & all_matrix[,4]);
  nACE=sum(only3 & all_matrix[,1] & all_matrix[,3] & all_matrix[,5]);
  nADE=sum(only3 & all_matrix[,1] & all_matrix[,4] & all_matrix[,5]);
  nBCD=sum(only3 & all_matrix[,2] & all_matrix[,3] & all_matrix[,4]);
  nBCE=sum(only3 & all_matrix[,2] & all_matrix[,3] & all_matrix[,5]);
  nBDE=sum(only3 & all_matrix[,2] & all_matrix[,4] & all_matrix[,5]);
  nCDE=sum(only3 & all_matrix[,3] & all_matrix[,4] & all_matrix[,5]);
  
  #elements overlapped in four data set.
  
  only4= apply(all_matrix,1,sum) == 4;
  nABCD=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3] &
              all_matrix[,4]);
  nABCE=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3] &
              all_matrix[,5]);
  nABDE=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,4] &
              all_matrix[,5]);
  nACDE=sum(only4 & all_matrix[,1] & all_matrix[,3] & all_matrix[,4] &
              all_matrix[,5]);
  nBCDE=sum(only4 & all_matrix[,2] & all_matrix[,3] & all_matrix[,4] &
              all_matrix[,5]);
  
  #elements overlapped in five data set.
  
  all5= apply(all_matrix,1,sum) == 5;
  nABCDE=sum(all5);
  
  #make the plot.
  
  elps=cbind(150*cos(seq(0,2*pi,len=1000)), 60*sin(seq(0,2*pi,len=1000)));
  
  relocate_elp=function(e, alpha, x, y){
    phi=(alpha/180)*pi;
    xr=e[,1]*cos(phi)+e[,2]*sin(phi);
    yr=-e[,1]*sin(phi)+e[,2]*cos(phi);
    xr=x+xr;
    yr=y+yr;
    return(cbind(xr, yr));
  }
  
  par(mar=c(1,1,1,1)); 
  plot(c(0, 400), c(0, 400), type="n", axes=F, ylab="", xlab="");
  
  polygon(relocate_elp(elps, 90,200, 250));
  polygon(relocate_elp(elps, 162,250, 220));
  polygon(relocate_elp(elps, 234,250, 150));
  polygon(relocate_elp(elps, 306,180, 125));
  polygon(relocate_elp(elps, 378,145, 200));
  
  #label the data set name.
  
  text(50, 280, colnames(all_matrix)[1]);
  text(170,400, colnames(all_matrix)[2]);
  text(350,300, colnames(all_matrix)[3]);
  text(350,20, colnames(all_matrix)[4]);
  text(55,10, colnames(all_matrix)[5]);
  
  #label the numbers
  
  text(61, 228, nA);
  text(194, 329, nB);
  text(321, 245, nC);
  text(290, 81, nD);
  text(132, 69, nE);
  
  text(146, 250, nAB, cex=0.7); 
  text(123, 188, nAC, cex=0.7); 
  text(275, 152, nAD, cex=0.7); 
  text(137, 146, nAE, cex=0.7); 
  text(243, 268, nBC, cex=0.7); 
  text(175, 267, nBD, cex=0.7); 
  text(187, 117, nBE, cex=0.7); 
  text(286, 188, nCD, cex=0.7); 
  text(267, 235, nCE, cex=0.7); 
  text(228, 105, nDE, cex=0.7); 
  
  text(148, 210, nABC,cex=0.7);
  text(159, 253, nABD,cex=0.7); 
  text(171, 141, nABE,cex=0.7); 
  text(281, 175, nACD,cex=0.7); 
  text(143, 163, nACE,cex=0.7); 
  text(252, 145, nADE,cex=0.7); 
  text(205, 255, nBCD,cex=0.7); 
  text(254, 243, nBCE,cex=0.7); 
  text(211, 118, nBDE,cex=0.7); 
  text(267, 211, nCDE,cex=0.7); 
  
  text(170, 231,nABCD,cex=0.7); 
  text(158, 169,nABCE,cex=0.7); 
  text(212, 139,nABDE,cex=0.7);
  text(263, 180,nACDE,cex=0.7); 
  text(239, 232,nBCDE,cex=0.7);
  
  text(204,190,nABCDE); 
}

