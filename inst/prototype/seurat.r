# https://nbisweden.github.io/single-cell_sib_scilifelab/session-differential-expression/differential_expression.html




FindAllMarkers <- function(object) {
    # marker_genes <- scran::findMarkers(
    #    object, 
    #    groups = colLabels(object), 
    #    test.type = "binom",
    #    lfc = 1, 
    #    direction = 'up')
    
    pwtt <- scran::pairwiseTTests(
            logcounts(object), 
            groups = colLabels(sce), 
            direction = "up")

    cbm_all <- scran::combineMarkers(
            de.lists = pwtt$statistics, pairs = pwtt$pairs,
            pval.type = "all"
        )


    lapply(names(cbm_all), function(cluster) {
            m <- cbm_all[[cluster]] |> as.data.frame()
            m$cluster <- cluster
            m$gene <- rownames(m)
            res <- m[, c("p.value", "FDR", "summary.logFC", "cluster", "gene")]
            return(res)
        }
    ) |> rbindlist()
}

sce$cluster_louvain_k10 <- colLabels(sce)

cowplot::plot_grid(scater::plotUMAP(sce, colour_by = "CD79A"),
                   scater::plotUMAP(sce, colour_by = "MS4A1"))

pp <- scater::plotHeatmap(sce, features = top10$gene,
                    columns = colnames(sce)[order(colLabels(sce))],
                    colour_columns_by = "cluster_louvain_k10", 
                    cluster_cols = F,
                    show_colnames = FALSE, cluster_rows = FALSE)









new.trend <- scran::modelGeneVarByPoisson(x = sce)
sce <- scran::denoisePCA(sce, technical = new.trend)


 

pbmc.markers <- FindAllMarkers(sce)


scater::plotExpression(sce, features = c("CD79A", "MS4A1"), 
                       x = "cluster_louvain_k10")


head(pbmc.markers)
sc_violin(sce, 'MS4A1', slot='logcounts') + ylim(0, 6)#+ scale_y_log10()
sc_dim(sce, reduction="UMAP")

sc_feature(sce, "CD3E", reduction='UMAP')
sc_feature(sce, "CD14", reduction='UMAP')

require(dplyr)
# pbmc.markers[order(pbmc.markers$FDR), ] %>%
mm |>
    group_by(cluster) %>%
    dplyr::filter(summary.logFC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

head(top10)

x <- assay(sce, 'scaled')[top10$gene,]
id <- Idents(sce)
colnames(x) <- names(id)


y <- mat2df(x) 
y$gene <- rownames(x)[y$row]
y$cell <- colnames(x)[y$col]
y$cluster <- id[y$cell]
y <- y[order(y$cluster), ]

gene_cls <- hclust(dist((x)))
cell_cls <- hclust(dist(t((x)))) 
gene_cls

y$gene <- factor(y$gene, levels = unique(gene_cls$labels[gene_cls$order]))
levels <- unique(cell_cls$labels[cell_cls$order])
y$cell <- factor(y$cell, levels = levels[order(id[levels])])
y$cell <- factor(y$cell, levels = names(sort(id)))
head(y)
library(ggfun)
p <- ggplot(y, aes(cell, gene, fill = value)) + 
    geom_tile()+ scale_fill_viridis_c(name = "Gene Expression") +
    theme_noxaxis() +
    scale_y_discrete(position = "right") +
    xlab(NULL) +
    ylab(NULL)
ggsave(p, file = 'xx.png')
p

library(Seurat)
FindMarkers.SingleCellExperiment <- function(object, slot, ...) {
    Seurat:::FindMarkers.Assay(object = object, slot = slot, ...)
}

GetAssayData.SingleCellExperiment <- function(object, slot) {
    assay(object, slot)
} 

xx =FindMarkers(sce, slot = 'logcounts', ident.1 = 5, ident.2 = 1)


pbmc.markers <- findMarkers(sce, pval.type = "all", lfc = 1, direction = "up")
head(pbmc.markers[[1]], 3)

cbm_all <- pbmc.markers
    mm =lapply(names(cbm_all), function(cluster) {
            m <- cbm_all[[cluster]] |> as.data.frame()
            m$cluster <- cluster
            m$gene <- rownames(m)
            res <- m[, c("p.value", "FDR", "summary.logFC", "cluster", "gene")]
            return(res)
        }
    ) |> rbindlist()


head(mm)
