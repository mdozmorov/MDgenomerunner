#' A function to convert V1 to V2 GenomeRunner results
#' 
#' Converts a log file of the results obtained with GenomeRunner V1 (standalone)
#' to the "matrix_PVAL/OR.txt" files compatible with GenomeRunner V2 (web)
#'
#' @param infile path to the log file. Required
#' @return generates two files, "matrix_PVAL.txt" and "matrix_OR.txt" in the same
#' filder as the log file
#' @export
#' @examples
#' \dontrun{
#' gr_getV1OddsRatioPvalMatrix("TFBS100_LOG.gr")
#' }
#' @note Courtesy to \link[https://github.com/AmyOlex]{Amy Olex}

gr_getV1OddsRatioPvalMatrix <- function(infile){
  
  ## First we have to do some raw file formatting
  ## Create the command to insert a header line for the read.delim function
  modfile <- paste(infile, "_modified", sep="")
  cmd <- paste("echo -e \"foi_name\tObserved\tExpected\tDiff\tp-val\tPCC\tObs/Tot\" | cat - ",infile," > ",modfile, sep="\t")
  ## execute the command and then remove the extraneous "-e " in the file (not sure why it puts it there, doing it on terminal doesn't have this issue)
  system(cmd)
  system(paste("cat ",modfile," | sed -e 's/^-e //' > tmp.txt && mv tmp.txt ",modfile, sep=""))
  
  # Now, load file after manually copying the header line into the file:
  v1 <- read.delim(modfile)
  
  # Get a list of GFs:
  GF <- subset(v1,grepl("^Features analyzed:.*",foi_name))
  
  #parse out GF name and total number
  parsed_GF <- as.data.frame(str_match(GF$foi_name, "Features analyzed: (.*) \\(Total (.*)\\)"))[,2:3]
  names(parsed_GF) <- c("GF", "tot")
  
  #add in column for line number
  parsed_GF$line <- as.numeric(row.names(GF))
  
  #add data start line number
  parsed_GF$dstart <- parsed_GF$line+9
  
  if(dim(GF)[1] == 1){
    parsed_GF$dend <- dim(v1)[1]
  } else { 
    #get data end numbers (next line - 2)
    tmp <- parsed_GF$line[2:length(parsed_GF$line)]-2
    parsed_GF$dend <- c(tmp, dim(v1)[1])
  }
  
  #get patients/disease features
  pdf <- as.character(v1[parsed_GF$dstart[1]:parsed_GF$dend[1],]$foi_name)
  pdf <- sub(".bed", "", pdf) # Remove ".bed", if any
  for (i in 1:length(pdf)) {
    pdf[i] <- tail(unlist(strsplit(pdf[i], "\\\\")), n=1) # Remove full path
  }
  
  data <- list()
  
  #add data to the list and calculate odds ratio and insert column
  for(i in 1:dim(parsed_GF)[1]){
    data[[i]] <- v1[parsed_GF$dstart[i]:parsed_GF$dend[i],]
    
    tmp <- data.frame(a=as.numeric(paste(data[[i]][,"Observed"])), c=as.numeric(paste(data[[i]][,"Expected"])) )
    tmp$b = as.numeric(paste(parsed_GF$tot[i])) - tmp$a
    tmp$d = as.numeric(paste(parsed_GF$tot[i])) - tmp$c
    ## Perform Fisher's exact test and store all the results
    tmp1 <- apply(tmp, 1, function(x) fisher.test(matrix(unlist(x), 2, 2)))
    # If OR confidence interval includes 1, then OR is not significant, set to 1
    data[[i]]$newOR <- sapply(tmp1, function(x) {ifelse(x$conf.int[1] < 1 & x$conf.int[2] > 1, 1, ifelse(x$estimate < 1, x$conf.int[2], x$conf.int[1]))})
    data[[i]]$newOR[ is.infinite(data[[i]]$newOR) ] <- .Machine$integer.max # Set infinite ORs to maximum number
    data[[i]]$newOR[ data[[i]]$newOR == 0 ] <- .Machine$double.eps # Set zero ORs to minimum number
    
    # Now modify the pvalues based on this new odds ratio value.
    data[[i]]$newPval <- sapply(tmp1, "[[", "p.value") # Store p-values
    data[[i]]$newPval[ data[[i]]$newPval < 1e-300 ] <- 1e-300 # Make 0 p-values to have some value, so they can be modified by OR
    data[[i]]$newPval[which(data[[i]]$newOR < 1)] = -1*data[[i]]$newPval[which(data[[i]]$newOR < 1)]
  }
  
  # Now we have all info we need to extract columns into matricies!
  # get the odds ratio data
  oddsRatios <- as.data.frame(lapply(data, "[[", "newOR"))
  row.names(oddsRatios) = pdf
  colnames(oddsRatios) = parsed_GF$GF
  # write odds ratio data to file
  ## Now save odds matrix to file
  write.table(as.data.frame((oddsRatios)), file=paste(dirname(infile), "matrix_OR.txt", sep="/"), quote=FALSE, sep="\t")
  
  # Now get the pvalues.  They are modified to have a negative sign IF the odds ratio is less than 1.  
  # If odds ratio is NaN or Inf or -Inf or NA then nothing is done to the pvalues.  
  # this happens when the contingency table data contains zeros.
  modPvals <- lapply(data, "[[", "newPval")
  modPvals <- lapply(modPvals, "paste")
  modPvals <- as.matrix(as.data.frame(lapply(modPvals, "as.numeric")))
  row.names(modPvals) = pdf
  colnames(modPvals) = parsed_GF$GF
  # Now save the modified pvalues to a file
  write.table(as.data.frame((modPvals)), file=paste(dirname(infile), "matrix_PVAL.txt", sep="/"), quote=FALSE, sep="\t")
}