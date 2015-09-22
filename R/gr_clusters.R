#' Define clusters
#' 
#' A function to define clusters by cutting a dendrogram by specified height. 
#' 
#' @param d a dendrogram object. Required.
#' @param height a height threshold to cut the dendrogram. Should be defined empirically. Required, default - 10.
#' @param minmembers a minimum number of samples to be considered as a cluster. Required, default - 3.
#' @param fileName save the results into a fileName. Default - none
#'
#' @return a list with 'eset.labels' and 'eset.groups' slots defining clustered labels/groups 
#' from the dendrogram cutting save the clustering order in a file
#' @export
#' @examples
#' \dontrun{
#' mtx.clust <- gr_clusters(d=h$Colv, height=1.5, minmembers=3, fileName="clustering.txt")
#' }
##
gr_clusters <- function(d, height = 10, minmembers = 3, fileName = NULL) {
  c <- cut(d, h = height)
  # Check the number of clusters, and the number of members.
  for (i in 1:length(c$lower)) {
    cat(paste("Cluster", formatC(i, width = 2, flag = "0"), sep = ""), "has", 
        formatC(attr(c$lower[[i]], "members"), width = 3), "members", "\n")
    if (!is.null(fileName)) {
      write.table(paste(i, t(labels(c$lower[[i]])), sep = "\t"), fileName, 
                  sep = "\t", quote = F, col.names = F, row.names = F, append = T)
    }
  }
  # Define Groups
  eset.labels <- character()  # Empty vector to hold cluster fileNames
  eset.groups <- numeric()  # Empty vector to hold cluster groups
  for (i in 1:length(c$lower)) {
    # Go through each cluster If the number of members is more than a minimum
    if (attr(c$lower[[i]], "members") > minmembers) {
      eset.labels <- append(eset.labels, labels(c$lower[[i]]))
      eset.groups <- append(eset.groups, rep(i, length(labels(c$lower[[i]]))))
    }
  }
  return(list(eset.labels = eset.labels, eset.groups = eset.groups))
}