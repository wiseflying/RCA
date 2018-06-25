#' Generate cell clusters
#' 
#' Hierarchical clustering and dynamic tree cutting. 
#' 
#' @param obj_in data object.
#' @param method can only be "hclust" (default). 
#' @param deepSplit_wgcna integer value indicating how deep the dendrogram should be cut. Values can be 0, 1 (default), 2, 3 and 4. 
#' @param min_group_Size_wgcna integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return data object.
#' @export
#' @examples
#' 
#' data_obj = cellClust(data_obj);
#' 
cellClust <- function(obj_in,method="hclust",deepSplit_wgcna=1,min_group_Size_wgcna=5)
{
  #### 1: reading input ####  
  fpkm_temp = obj_in$fpkm_for_clust;
  
  if (!require(flashClust)) install.packages("flashClust",repos = "http://cran.us.r-project.org") 
  require(flashClust)
  #if (!require(WGCNA)){
  #  source("http://bioconductor.org/biocLite.R");
  #  biocLite(c("impute", "GO.db", "preprocessCore")); 
  #  install.packages("WGCNA");
  #}
  #require(WGCNA)
  require(dynamicTreeCut)
  labels_to_colors <-function (labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey", 
    commonColorCode = TRUE) 
{
    if (is.null(colorSeq)) 
        colorSeq = standardColors()
    if (is.numeric(labels)) {
        if (zeroIsGrey) 
            minLabel = 0
        else minLabel = 1
        if (any(labels < 0, na.rm = TRUE)) 
            minLabel = min(c(labels), na.rm = TRUE)
        nLabels = labels
    }
    else {
        if (commonColorCode) {
            factors = factor(c(as.matrix(as.data.frame(labels))))
            nLabels = as.numeric(factors)
            dim(nLabels) = dim(labels)
        }
        else {
            labels = as.matrix(as.data.frame(labels))
            factors = list()
            for (c in 1:ncol(labels)) factors[[c]] = factor(labels[, 
                c])
            nLabels = sapply(factors, as.numeric)
        }
    }
    if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
        nRepeats = as.integer((max(labels) - 1)/length(colorSeq)) + 
            1
        warning(paste("labels2colors: Number of labels exceeds number of avilable colors.", 
            "Some colors will be repeated", nRepeats, "times."))
        extColorSeq = colorSeq
        for (rep in 1:nRepeats) extColorSeq = c(extColorSeq, 
            paste(colorSeq, ".", rep, sep = ""))
    }
    else {
        nRepeats = 1
        extColorSeq = colorSeq
    }
    colors = rep("grey", length(nLabels))
    fin = !is.na(nLabels)
    colors[!fin] = naColor
    finLabels = nLabels[fin]
    colors[fin][finLabels != 0] = extColorSeq[finLabels[finLabels != 
        0]]
    if (!is.null(dim(labels))) 
        dim(colors) = dim(labels)
    colors
}
  
  
  
  
  
  
  
  #### 2: choosing method ####
    if (method == "hclust"){
        d = as.dist(1-cor(fpkm_temp,method="pearson"));
        cellTree = flashClust(d,method = "average");     
        dynamicGroups = cutreeDynamic(dendro = cellTree,distM = as.matrix(d),deepSplit = deepSplit_wgcna,
                                      pamStage = FALSE,minClusterSize= min_group_Size_wgcna);
        dynamicColors = labels2colors(dynamicGroups);
        group_labels = matrix(dynamicGroups,nrow = length(names(fpkm_temp)),ncol = 1,list(names(fpkm_temp),c("groupLabel")),byrow=FALSE)
        group_labels = as.data.frame(group_labels)
        group_labels_color = cbind(group_labels,dynamicColors);
    }
  #### 3: writing output ####    
  obj_out = append(obj_in,
              list("d" = d,
              "cellTree" = cellTree,
              "dynamicColors" = dynamicColors,
              "group_labels_color" = group_labels_color
              )
  )
  
  return(obj_out)
  
}
