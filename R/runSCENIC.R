#' wrapper function for running SCENIC on default parameter. 
#' 
#' https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
#' 
#' @param sce A singlecellExperiment object.
#' @param species The species of the data, one of human, mouse or drosophila_melanogaster.
#' @param nCores The number of cores to use for parallel processing.
#' @param dbDir The directory path to the cisTarget databases files.
#' @param assay_name the name of assay.
#' @param minCountsPerGene The minimum number of counts per gene required for inclusion in the analysis.
#' @param minSamples The minimum proportion of samples in which a gene must be detected to be included in the analysis.
#' 
#' @return A matrix of regulon activity scores for each cell.
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @export 
runSCENIC <- function(sce,species,nCores,dbDir,
                      assay_name = "logcounts",
                      minCountsPerGene = 3,
                      minSamples = 0.01){
  dir.create("int")
  # 设置物种和TF数据库路径
  org <- switch(species, 
                "human" = "hgnc", 
                "mouse" = "mgi", 
                "drosophila_melanogaster" = "dmel", 
                stop("Unsupported species"))
  dbs <- SCENIC::defaultDbNames[[org]]
  
  if(! "logcounts" %in% SummarizedExperiment::assayNames(sce)){
    sce <- Seurat::NormalizeData(sce)
  }
  
  # 构建scenicOptions对象
  scenicOptions <- SCENIC::initializeScenic(org = org, dbDir = dbDir, dbs = dbs, nCores = nCores)
  scenicOptions@inputDatasetInfo$cellInfo <- SummarizedExperiment::colData(sce)
  
  # 过滤
  exprMat = SummarizedExperiment::assay(sce, assay_name)
  genesKept <- SCENIC::geneFiltering(exprMat, scenicOptions = scenicOptions,
                                     minCountsPerGene = minCountsPerGene * ncol(sce),
                                     minSamples = minSamples * ncol(sce)) 
  exprMat_filtered <- exprMat[genesKept, ]
  # 计算相关性
  SCENIC::runCorrelation(exprMat_filtered, scenicOptions)
  
  # 运行GENIE3
  set.seed(123)
  SCENIC::runGenie3(exprMat_filtered, scenicOptions)
  
  # 构建GRN和打分
  scenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- SCENIC::runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMat)
 
  return(scenicOptions)
}

