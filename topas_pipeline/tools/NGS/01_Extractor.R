library('MultiAssayExperiment')
library("RaggedExperiment")
library('GenomicRanges')
library(data.table)
library('magrittr')
library(dplyr)



#' Extract Data Across Modalities by Genomic Overlap
#'
#' This function retrieves overlapping genomic features from a `MultiAssayExperiment`-like object
#' for specified data modalities (e.g., SNV, INDEL, CNV, FUSION), based on input query regions (e.g., gene coordinates).
#'
#' @param dataObj A list-like object (e.g., MultiAssayExperiment or named list) containing genomic datasets,
#'   where each element (e.g., "snv", "cnv", etc.) is coercible to a `GRangesList`.
#' @param query A `GRanges` object representing the genomic regions of interest (default is genes from `gencode19_gns_lite`
#'   filtered by names in `all_genes`).
#' @param dts A character vector of data types/modalities to extract from `dataObj`. Each name should correspond
#'   to an element in `dataObj`. Default is `c("snv", "indel", "cnv", "fusion")`.
#'
#' @return A named list with one element per modality in `dts`. Each element is itself a named list where
#'   names are gene names and values are `GRanges` objects representing overlapping features.
#'
#' @examples
#' # Example: Extract overlapping features for selected genes
#' result <- get_data_for_all_modalities(dataObj = mae_object,
#'                                       query = some_gene_regions,
#'                                       dts = c("snv", "cnv"))
#'
#' @import GenomicRanges
#' @export
get_data_for_all_modalities <- function(dataObj, 
                                        query = gencode19_gns_lite[gencode19_gns_lite$gene_name %in% all_genes],
                                        dts = c("snv", "indel", "cnv", "fusion")
                                        ){
  
   tmp = setNames(lapply(dts, function(dt) {
     # convert MAE type dt to GRL
     obj_grl = as(dataObj[[dt]], "GRangesList")
     res = list()
     # iterate over queries
     for(i in 1:length(query)) {
       obj_gr = unlist(obj_grl)
       res[[query$gene_name[i]]] = subsetByOverlaps(obj_gr, query[i])
     }
     res
   }), nm=dts)
  
  return(tmp)
}




#' Extract Per-Protein Data for a Given Modality
#'
#' This function extracts patient-level data for each protein (gene) from a specific modality 
#' previously returned by `get_data_for_all_modalities()`. It returns a combined data frame 
#' with information from all requested proteins.
#'
#' @param tmp A named list produced by `get_data_for_all_modalities()`, where each element corresponds 
#'   to a modality (e.g., "snv", "cnv") and contains a list of `GRanges` objects per gene.
#' @param modality_name A character string specifying the modality name to extract from `tmp`
#'   (e.g., "snv", "cnv", "fusion").
#' @param list_genes A character vector of gene (or protein) names to extract from the selected modality.
#'
#' @return A `data.table::data.table` (or `data.frame` if `data.table` is not loaded) with rows corresponding 
#'   to genomic features per patient and per protein. Columns include:
#'   - `patient`: the patient or sample name
#'   - All metadata columns from the corresponding `GRanges` object
#'   - `protein_name`: the name of the gene/protein queried
#'
#' @details
#' The function safely handles missing or malformed data using `tryCatch`, printing warnings or errors 
#' for specific genes but continuing execution.
#'
#' @examples
#' # Assuming 'tmp' is the result of get_data_for_all_modalities(...)
#' gene_list <- c("TP53", "EGFR")
#' snv_data <- get_data_per_protein_per_modality(tmp, modality_name = "snv", list_genes = gene_list)
#'
#' @import data.table
#' @export
get_data_per_protein_per_modality <- function(tmp, modality_name, list_genes) {
  df <- tmp[[modality_name]]
  
  get_data_per_protein <- function(df, protein_name) {
    tryCatch({
      df_protein <- data.frame(
        patient = names(df[[protein_name]]),
        df[[protein_name]]
      )
      df_protein$protein_name <- protein_name
      return(df_protein)
    }, error = function(e) {
      message("Error processing protein: ", protein_name)
      message("Details: ", e$message)
      return(data.frame())
    }, warning = function(w) {
      message("Warning processing protein: ", protein_name)
      message("Details: ", w$message)
      return(data.frame())
    })
  }
  list_dfs = lapply(list_genes, function(x) get_data_per_protein(df, x))
  return(rbindlist(list_dfs))
}


### Running the script

fp_csv = read.csv('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2025.07.22_CJ_paper_cohort/annot_fp.csv')
all_genes = fp_csv$Gene.names %>% strsplit(";") %>% unlist() %>% unique()

load('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/gencode19_gns_lite.RData')
load('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/dataCHORDOMA_240515.RData')
load('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/dataMASTER_250526215122.RData')


tmp = get_data_for_all_modalities(dataMASTER)
saveRDS(tmp,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/master.rds' )
tmp <- readRDS('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/master.rds')
fusion = get_data_per_protein_per_modality(tmp,'fusion',all_genes)
write.csv(fusion,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/fusion_df.csv')
cnv_df = get_data_per_protein_per_modality(tmp,'cnv',all_genes)
write.csv(cnv_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/cnv_df.csv')
snv_df = get_data_per_protein_per_modality(tmp,'snv',all_genes)
write.csv(snv_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/snv_df.csv')
indel_df = get_data_per_protein_per_modality(tmp,'indel',all_genes)
write.csv(indel_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome/indel_df.csv')


tmp = get_data_for_all_modalities(dataCHORDOMA)
saveRDS(tmp,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/chordoma.rds' )
tmp <- readRDS('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/chordoma.rds')
fusion = get_data_per_protein_per_modality(tmp,'fusion',all_genes)
write.csv(fusion,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/fusion_df.csv')
cnv_df = get_data_per_protein_per_modality(tmp,'cnv',all_genes)
write.csv(cnv_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/cnv_df.csv')
snv_df = get_data_per_protein_per_modality(tmp,'snv',all_genes)
write.csv(snv_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/snv_df.csv')
indel_df = get_data_per_protein_per_modality(tmp,'indel',all_genes)
write.csv(indel_df,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM/indel_df.csv')
