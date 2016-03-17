#' @title Extract FastQC summary tables
#' 
#' @description \code{fastqc} function parses all the various summary statistics produced by \code{FastQC} and returns a nested list data structure corresponding to each sample.
#' 
#' @param sample_info Full path to file containing details about \code{samples} and their \code{paths}. See \code{Details} section.
#' 
#' @details The file provided to \code{sample_info} argument should contain at least these three columns: 
#' \itemize{
#'   \item \code{sample} - contains the \code{sample} name.
#'   \item \code{pair} - in case of paired end reads, \code{1} or \code{2} corresponding to first and second pair, and in case of single end reads, \code{NA}.
#'   \item \code{path} - full path to the fastqc summary report (\code{.txt} file) for each sample.
#' }
#' 
#' It can also optionally contain a \code{group} column. If present, then the plot generated with fastqc would be colored / facetted according to samples belonging to the same group.
#' @seealso \code{\link{plot_dup_stats}}, \code{\link{plot_gc_stats}}, \code{\link{plot_sequence_quality}}, \code{\link{plot_total_sequence_stats}}
#' @export
#' @examples
#' \dontrun{
#'   require(ggfastqc)
#'   obj = fastqc("~/Documents/oa/omics/CRC_AZA/fastqc-sample-file.tsv")
#' }
fastqc <- function(sample_info) {
    info = fread(sample_info)
    cols = c("sample", "pair", "path")
    rest = setdiff(cols, names(info))
    if (length(rest)) stop("Columns [", paste(rest, collapse=","), "] not found in file provided to argument 'sample_info'.")
    pair_uniq = sort(unique(info[["pair"]]))
    if (length(pair_uniq) == 1L && !is.na(pair_uniq)) {
        if (pair_uniq != 1L) stop("Pair column must be NA for single end reads, and 1 or 2 for paired end reds.")
        else warning("Pair column must be NA for single end reads.")
    } else if ( (length(pair_uniq) == 2L && !identical(pair_uniq, 1:2)) || !length(pair_uniq) %in% 1:2 ) {
        stop("Pair column must be NA for single end reads, and 1 or 2 for paired end reads.")
    }
    if (is.null(info[["group"]])) info[, group := ""]
    samples = if (length(pair_uniq) == 1L) info[["sample"]] else info[, paste(sample, pair, sep="_")]
    ans = info[, .(tables = .(extract_summary_tables(readLines(path), sample, pair, group))), by=1:nrow(info)][["tables"]]
    setattr(ans, 'class', 'fastqc')
    setattr(ans, 'names', samples)
}

#' @title Extract all summary tables for a particular sample
#' 
#' @description This is the workhorse that powers \code{fastqc}. Given a \code{sample} and \code{path}, \code{extract_summary_tables} extracts all summary statistics and returns them as a list of \code{data.table}s.
#' 
#' @param doc The file corresponding to \emph{that \code{sample}} read in as such (using \code{readLines}).
#' @param sample Sample name
#' @param pair  Read pair
#' @param group     Group / condition this sample belongs to
extract_summary_tables <- function(doc, sample, pair, group) {

    extract <- function(i, j) {
        if (j-i-2L <= 0L) {
             message("Table at start index ", i, " seems to contain only one (header) row, therefore returning empty data.table.")
             return(setDT(list())[]) # data.table:::null.data.table()
        }
        hdr_row = unlist(strsplit(doc[i+1L], "\t", fixed=TRUE))
        bdy_row = strsplit(doc[(i+2L):(j-1L)], "\t", fixed=TRUE)
        ans = setDT(transpose(bdy_row))
        if (ncol(ans) != length(hdr_row)) 
            hdr_row = c(hdr_row, paste("V", seq_len(ncol(ans)-length(hdr_row)), sep=""))
        setnames(ans, hdr_row)
        ans[, c("sample_group", "pair") := .(sample, pair)]
        if (!is.na(pair)) ans[, sample_name := paste(sample, pair, sep="_")]
        ans[, group := group]
        movecols = c("group", "sample_group", "sample_name", "pair")
        setcolorder( ans, c( movecols, head(names(ans), -length(movecols)) ) )
    }

    idx  = grep("^>>", doc)
    idx1 = idx[c(TRUE, FALSE)]
    idx2 = idx[c(FALSE, TRUE)]
    hdr  = gsub("^>>(.*)\t.*$", "\\1", doc[idx1])

    tables = lapply(seq_along(idx1), function(i) extract(idx1[i], idx2[i]))
    setattr(tables, 'names', hdr)
}

# deprecated, but kept for now
split_sample <- function(sample_name) {
    splits = tstrsplit(sample_name, "_(?=[12]$)", perl=TRUE)
    if (length(splits) == 1L)
        splits = c(splits, list(1L))
    if (length(splits) != 2L)
        stop("'sample_name' argument should be of the form <samplename>_<pair> for paired end and <samplename> for single end Fastq files, with no '_' elsewhere.")
    splits[[2]] = as.integer(splits[[2]])
    splits
    # x[, c("sample_name", "pair") := splits]
    # setcolorder(x, c("sample_name", "pair", head(names(x), -2L)))
}
