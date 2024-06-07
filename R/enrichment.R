## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: gene set enrichment analyses

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(ggplot2)

#' Gene set enrichment analysis on Gene Ontology terms
#'
#' A simple wrapper around [clusterProfiler::gseGO()].
#'
#' @param scores Named sorted list of scores for genes
#' @param ontology Subontology to use, or "ALL"
#' @param db Organism database
#' @param pvalue.cutoff p-value cut-off (after multiple testing adjustment)
#' @param key.type Key type used to map names of `scores` to entries in `db`
#' @param ... Further parameters passed to [gseGO()]
#' @return GSEA results
gsea.go <- function(scores, ontology=c("BP", "CC", "MF", "ALL"), db=org.Hs.eg.db,
                    pvalue.cutoff=0.05, key.type="SYMBOL", ...) {
  gseGO(scores, ontology[1], db, keyType=key.type, eps=0, pvalueCutoff=pvalue.cutoff, ...)
}


#' Gene set enrichment analysis on MSigDB gene sets
#'
#' @param scores Named sorted list of scores for genes
#' @param gene.sets Data frame of gene sets (as returned by [msigdbr()])
#' @param pvalue.cutoff p-value cut-off (after multiple testing adjustment)
#' @param ... Further parameters passed to [GSEA()]
#' @return GSEA results
gsea.msigdb <- function(scores, gene.sets=msigdbr(category="H"), pvalue.cutoff=0.05, ...) {
  GSEA(scores, TERM2GENE=gene.sets[, c("gs_name", "human_gene_symbol")],
       TERM2NAME=unique(gene.sets[, c("gs_name", "gs_description")]), pvalueCutoff=pvalue.cutoff, eps=0, ...)
}


#' Gene set enrichment analysis on Reactome pathways
#'
#' The underlying function [ReactomePA::gsePathway()] expects genes to be identified by Entrez IDs.
#' Hence a mapping between the names in `scores` (presumed gene symbols) and Entrez IDs is passed via `mapping`.
#' Entrez IDs of core genes in the results (component `@result`) are automatically replaced if `replace.ids` is `TRUE`.
#'
#' @param scores Named sorted list of scores for genes
#' @param mapping Data frame with mapping between gene symbols and Entrez IDs (as returned by [bitr()])
#' @param replace.ids Replace Entrez IDs in results with gene symbols?
#' @param pvalue.cutoff p-value cut-off (after multiple testing adjustment)
#' @param ... Further parameters passed to [gsePathway()]
#' @return GSEA results
gsea.reactome <- function(scores, mapping=bitr(names(scores), "SYMBOL", "ENTREZID", org.Hs.eg.db, FALSE),
                          replace.ids=TRUE, pvalue.cutoff=0.05, ...) {
  if (!is.null(mapping))
    names(scores) <- mapping[match(names(scores), mapping$SYMBOL), "ENTREZID"]
  gse.res <- gsePathway(scores, eps=0, pvalueCutoff=pvalue.cutoff, ...)
  if (replace.ids && (nrow(gse.res) > 0)) {
    ids <- strsplit(gse.res@result$core_enrichment, "/")
    genes <- lapply(ids, function(x) mapping[match(x, mapping$ENTREZID), "SYMBOL"])
    gse.res@result$core_enrichment <- sapply(genes, paste, collapse="/")
  }
  gse.res
}


#' Create a dot plot of enrichment results
#'
#' `(+)` or `(-)` are appended to term descriptions to indicate enrichment among the top or bottom end of the gene list, respectively.
#'
#' @param gse.res Gene set enrichment result (`gseaResult` object)
#' @param n Number of hits to show
#' @param label_format Space allocated to gene set labels
#' @param ... Further arguments passed to [dotplot()]
#' @return Dot plot
dotplot.direction <- function(gse.res, n=15, label_format=40, ...) {
  if (nrow(gse.res) == 0) {
    warning("No GSEA results to plot")
    return(ggplot())
  }
  updown <- ifelse(gse.res$enrichmentScore < 0, "(-)", "(+)")
  gse.res@result$Description <- paste(sub("\\.$", "", gse.res@result$Description), updown)
  dotplot(gse.res, order="pvalue", decreasing=FALSE, showCategory=n, label_format=label_format, ...)
}


#' Run gene set enrichment analyses for different sources of gene sets
#'
#' Performs GSEA based on GO terms ("biological process" and "molecular function" ontologies), Reactome pathways, and MSigDB Hallmark gene sets.
#'
#' @param scores Named list of scores for genes
#' @param out.prefix Path and filename prefix for output files (extensions will be appended)
#' @param plot Generate PDF file with plots?
#' @param reactome.mapping Data frame with mapping between gene symbols and Entrez IDs (as returned by [bitr()])
#' @param ... Further parameters passed to the `gsea...()` functions
#' @return List of GSEA results
gsea.all <- function(scores, out.prefix="GSEA_results", plot=TRUE,
                     reactome.mapping=bitr(names(scores), "SYMBOL", "ENTREZID", org.Hs.eg.db, FALSE), ...) {
  scores <- na.omit(sort(scores, decreasing=TRUE))

  message("GSEA: Gene Ontology (Biological Process)...")
  gse.go.bp <- gsea.go(scores, "BP", ...)
  message("GSEA: Gene Ontology (Molecular Function)...")
  gse.go.mf <- gsea.go(scores, "MF", ...)
  message("GSEA: Reactome pathways...")
  gse.reactome <- gsea.reactome(scores, mapping=reactome.mapping, ...)
  message("GSEA: Hallmark gene sets...")
  gse.hallmark <- gsea.msigdb(scores, ...)

  if (nchar(out.prefix) > 0) {
    write.csv(gse.go.bp@result, paste0(out.prefix, "_GO-BP.csv"))
    write.csv(gse.go.mf@result, paste0(out.prefix, "_GO-MF.csv"))
    write.csv(gse.reactome@result, paste0(out.prefix, "_Reactome.csv"))
    write.csv(gse.hallmark@result, paste0(out.prefix, "_Hallmark.csv"))

    if (plot) {
      pdf(paste0(out.prefix, ".pdf"))
      ## theme_update(plot.title=element_text(hjust=0.5, face="bold")) # not effective
      print(dotplot.direction(gse.go.bp, label_format=40) +
            labs(title="Gene Ontology (Biological Process)") +
            theme(plot.title=element_text(hjust=0.5, face="bold")))
      print(dotplot.direction(gse.go.mf, label_format=40) +
            labs(title="Gene Ontology (Molecular Function)") +
            theme(plot.title=element_text(hjust=0.5, face="bold")))
      print(dotplot.direction(gse.reactome, label_format=40) +
            labs(title="Reactome Pathways") +
            theme(plot.title=element_text(hjust=0.5, face="bold")))
      print(dotplot.direction(gse.hallmark, label_format=55) +
            labs(title="Hallmark Gene Sets") +
            theme(plot.title=element_text(hjust=0.5, face="bold")))
      dev.off()
    }
  }

  invisible(list(GO.BP=gse.go.bp, GO.MF=gse.go.mf, Reactome=gse.reactome, Hallmark=gse.hallmark))
}
