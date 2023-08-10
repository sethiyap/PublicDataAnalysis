#' get_fpkm
#'
#' @param gff
#' @param count_file
#' @param percent_mapped_file
#' @param percent_map_cutoff
#'
#' @return
#' @export
#'
#' @examples
#' gff_file <- "~/Documents/CDK7_project/RNASeq/AGRF_CAGRF220611046_HGHV5DRX2/GENOME/FungiDB-59_CneoformansH99.gff"
#'

get_fpkm <- function(gff, count_file, percent_mapped_file, percent_map_cutoff=70){

  gff <- GenomicFeatures::makeTxDbFromGFF(file = gff_file)

  # determine gene length to calculate fpkm
  gg <- GenomicFeatures::genes(gff)

  # tibble of gene length
  gene_length <- tidyr::as_tibble(lengths(gg), rownames = "genes") %>%
    dplyr::rename("length"="value")


  # read the countfile obtained from rsubread
  count_file <- readr::read_delim(count_file, delim="\t", col_names = TRUE) %>%
    dplyr::inner_join(gene_length,by=c("gene_name"="genes"))

  # join with mapping percentage
  mapping_percentage <- readr::read_delim(percent_mapped_file, delim="\t", col_names = TRUE) %>%
    dplyr::mutate(percent=gsub(pattern = "%", replacement = "", x = percent))

  # join and filter samples based on mapping percentage
  count_mat_tibble <- count_file %>%
    tidyr::gather(samples, value, -gene_name, -length) %>%
    tidyr::nest(cols=c(gene_name, value, length)) %>%
    dplyr::inner_join(mapping_percentage, by=c("samples"="names")) %>%
    dplyr::filter(percent > percent_map_cutoff)


  calculate_fpkm <- count_mat_tibble %>%
    dplyr::mutate(fpkm = purrr::map(cols, function(ii){



      rpm <- ii$value*1e9

      rpkm <- rpm/(ii$length*sum(ii$value))

    tibble::tibble(gene_name = ii$gene_name, fpkm = rpkm) })) %>%
    dplyr::select(c(samples,fpkm)) %>%
    tidyr::unnest(cols = fpkm) %>%
    tidyr::pivot_wider(names_from  = samples, values_from = fpkm)

  return(calculate_fpkm)

}



