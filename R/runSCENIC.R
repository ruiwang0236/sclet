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
#' @importFrom 
#' 
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
  data(defaultDbNames)
  dbs <- defaultDbNames[[org]]
  
  if(! "logcounts" %in% assayNames(sce)){
    sce <- NormalizeData(sce)
  }
  
  # 构建scenicOptions对象
  scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, nCores = nCores)
  scenicOptions@inputDatasetInfo$cellInfo <- colData(sce)
  
  # 过滤
  exprMat = assay(sce, assay_name)
  genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions,
                             minCountsPerGene = minCountsPerGene * ncol(sce),
                             minSamples = minSamples * ncol(sce)) 
  exprMat_filtered <- exprMat[genesKept, ]
  # 计算相关性
  runCorrelation(exprMat_filtered, scenicOptions)
  
  # 运行GENIE3
  set.seed(123)
  runGenie3(exprMat_filtered, scenicOptions)
  
  # 构建GRN和打分
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)
  
  # 读取第三步打分结果
  regulonAUC <- readRDS('int/3.4_regulonAUC.Rds')
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  
  # 获得regulon打分矩阵
  regAct <- getAUC(regulonAUC)
  
  return(regAct)
}

