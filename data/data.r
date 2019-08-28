#' Foliar endophytic fungi across the Hawaiian Archipelago
#'
#' A dataset containing an OTU table (species-by-site), environmental
#' metadata, and host plant phylogeny.
#'
#' @format A list containing 3 objects:
#' \item{zotutable}{
#'   data.frame object where each row is a sample and each column is a
#'   fungal zOTU. Rownames are sample IDs.
#' }
#' \item{metadata}{
#'   data.frame object containing environmental metadata for samples in
#'   zotutable. SampleID column of metadata matches rownames of zotutable.
#' }
#' \item{supertree}{
#'   Phylogenetic tree containing all host plant genera in PlantGenus
#'   column of metadata.
#' }
#' @source no source yet.
"endophyte"