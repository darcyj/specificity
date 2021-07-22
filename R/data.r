#' Foliar endophytic fungi across the Hawaiian Archipelago
#'
#' A dataset containing an OTU table (species-by-site), environmental
#' metadata, and host plant phylogeny.
#'
#' @format A list containing 4 objects:
#' \describe{
#'   \item{otutable:}{
#'     data.frame object where each row is a sample and each column is a
#'     fungal OTU (actually ASV from DADA2). Rownames are sample IDs.
#'   }
#'   \item{metadata:}{
#'     data.frame object containing environmental metadata for samples in
#'     otutable. SampleID column of metadata matches rownames of otutable.
#'   }
#'   \item{supertree:}{
#'     Phylogenetic tree containing all host plant genera in PlantGenus
#'     column of metadata.
#'   }
#'   \item{taxonomy:}{
#'     UNITE (v8) taxonomy for OTUs in otutable. Assigned using BLAST, e-value (eval) and UNITE best hit (hit) are given as well. 
#'   }
#' }
#' @source Darcy et al. (2020) Fungal communities living within leaves 
#'   of native Hawaiian dicots are structured by landscape‐scale 
#'   variables as well as by host plants. Mol Ecol 29:3102-3115 
#'   https://doi.org/10.1111/mec.15544
#' 
"endophyte"

