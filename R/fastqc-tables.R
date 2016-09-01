#' @title Extract FastQC summary tables
#' 
#' @description \code{fastqc} function parses all the summary statistics of 
#' reports produced by \code{FastQC} tool and returns a \code{data.table} with 
#' two columns: \code{param} and \code{value}. 
#' 
#' Each row of the \code{value} column contains the data corresponding to that 
#' \code{param}, and is itself 
#' a \code{data.table}.
#' 
#' @param sample_info Full path to file containing details about 
#' \code{samples} and their \code{paths}.
#' 
#' @details The file provided to \code{sample_info} argument should contain 
#' at least these three columns: 
#' \itemize{
#'   \item \code{sample} - contains the \code{sample} name.
#'   \item \code{pair} - in case of paired end reads, \code{1} or \code{2} 
#'   corresponding to first and second pair, and in case of single end reads, 
#'   \code{NA}.
#'   \item \code{path} - full path to the fastqc summary report (\code{.txt} 
#'   file) for each sample.
#'    
#'   If just the file name (\code{.txt}) is provided, it is assumed that the 
#'   file is in the same folder as the input file provided to 
#'   \code{sample_info} argument.
#' }
#' 
#' It can also optionally contain a \code{group} column. If present, the plots 
#' generated will take it into account and \code{color} / \code{facet} 
#' accordingly.
#' 
#' @seealso \code{\link{plot_dup_stats}} \code{\link{plot_gc_stats}} 
#' \code{\link{plot_sequence_quality}} \code{\link{plot_total_sequence_stats}}
#' @return An object of class \code{fastqc} which inherits from 
#' \code{"data.table"}, with two columns: \code{param} and \code{value}, where 
#' \code{value} is a list of \code{"data.table"}s.
#' @export
#' @examples
#' path = system.file("tests/fastqc-sample", package="ggfastqc")
#' obj = fastqc(sample_info = file.path(path, "annotation.txt"))
fastqc <- function(sample_info) {
    info = fread(sample_info, colClasses=list(character=c("sample")))
    cols = c("sample", "pair", "path")
    rest = setdiff(cols, names(info))
    if (length(rest)) stop("Columns [", paste(rest, collapse=","), 
        "] not found in file provided to argument 'sample_info'.")
    pair = group = path = NULL
    pair_uniq = sort(unique(info[["pair"]]))
    if (length(pair_uniq) == 1L && !is.na(pair_uniq)) {
        if (pair_uniq != 1L) stop("Pair column must be NA for single end ", 
                                "reads, and 1 or 2 for paired end reds.")
        else warning("Pair column must be NA for single end reads.")
    } else if ( (length(pair_uniq) == 2L && 
                !identical(pair_uniq, 1:2)) || 
                !length(pair_uniq) %in% 1:2 ) {
        stop("Pair column must be NA for single end reads, ", "
                and 1 or 2 for paired end reads.")
    }
    if (is.null(info[["group"]])) info[, "group" := ""]
    samples = if (length(pair_uniq) == 1L) info[["sample"]] 
                else info[, paste(sample, pair, sep="_")]
    # set paths correctly, if necessary
    idx = which(info[["path"]] == basename(info[["path"]]))
    if (length(idx)) 
        info[(idx), "path" := file.path(dirname(sample_info), path)][]
    ans = info[, list(tables = list(extract_summary_tables(readLines(path), 
                sample, pair, group))), by=1:nrow(info)][["tables"]]
    fans = lapply(seq_along(ans[[1L]]), function(i) 
            rbindlist(lapply(ans, `[[`, i)))
    fans = setDT(list(param = names(ans[[1L]]), value=fans))
    setattr(fans, 'class', c('fastqc', class(fans)))
}

#' @title Extract all summary tables for a particular sample
#' 
#' @description \bold{NOTE:} For internal use only.
#' 
#' This is the workhorse that powers \code{fastqc}. Given a \code{sample} and 
#' \code{path}, \code{extract_summary_tables} extracts all summary statistics 
#' and returns them as a list of \code{data.table}s.
#' 
#' @return A list of data.tables.
#' @param doc The file corresponding to \emph{that \code{sample}} read in 
#' as such (using \code{readLines}).
#' @param sample Sample name.
#' @param pair  Read pair.
#' @param group Group / condition this sample belongs to.
extract_summary_tables <- function(doc, sample, pair, group) {

    fix_style <- function(v) {
        v = gsub("[ ]*[/][ ]*", ".", tolower(v))
        v = gsub("[ ]+|[-]", "_", v)
        v = gsub("^[#]", "", v)
        v = gsub("^[%]", "percent_", v)
    }

    extract <- function(i, j) {
        stopifnot(length(i) == 1L, length(j) == 1L)
        # parse and create data.table
        hdr_row = unlist(strsplit(doc[i], "\t", fixed=TRUE))
        if (is.na(i) || is.na(j)) {
            return (setDT(list())) # null data.table so that rbindlist skips it
        } else if (j-i <= 0L) {
            # the entire chunk contains only this header
            # make it a 1-col data.table
            bdy_row = hdr_row
            ans = data.table(bdy_row[-1L])
            setnames(ans, bdy_row[1L])
        } else {
            bdy_row = strsplit(doc[(i+1L):(j)], "\t", fixed=TRUE)
            ans = setDT(transpose(bdy_row))
            if (ncol(ans) != length(hdr_row)) 
                hdr_row = c(hdr_row, paste("V", seq_len(ncol(ans) - 
                                    length(hdr_row)), sep=""))
            setnames(ans, hdr_row)
        }
        # make sure numeric/integer cols are of numeric/integer types
        ans[, names(ans) := lapply(.SD, type.convert, as.is=TRUE), 
                .SDcols=names(ans)]
        # fix style
        chr_cols = names(ans)[sapply(ans, is.character)]
        if (length(chr_cols)) {
            ans[, (chr_cols) := lapply(.SD, fix_style), .SDcols=chr_cols]
        }
        # add sample_group, pair and group cols
        ans[, c("sample_group", "pair") := list(sample, pair)]
        if (!is.na(pair)) ans[, "sample_name" := paste(sample, pair, sep="_")]
        ans[, "group" := group]
        # lowercase colnames and fix them to be proper
        setnames(ans, fix_style(names(ans)))
        # convert sequences to lowercase as well
        if ("sequence" %in% names(ans)) {
            ans[, "sequence" := tolower(sequence)][]
        }
        # reorder cols
        movecols = c("group", "sample_group", "sample_name", "pair")
        setcolorder( ans, c(movecols, head(names(ans), -length(movecols))) )
    }

    idx  = grep("^>>", doc)
    idx1 = idx[c(TRUE, FALSE)]
    idx2 = idx[c(FALSE, TRUE)]
    hdr  = gsub("^>>(.*)\t.*$", "\\1", doc[idx1])
    hdr  = gsub("[ ]+", "_", tolower(hdr))

    # getting around the messy fastqc data format..
    eq_len = (idx2-idx1 == 1L)
    idx1 = ifelse(eq_len, NA_integer_, idx1+1L)
    idx2 = ifelse(eq_len, NA_integer_, idx2-1L)

    dup_hdr = "sequence_duplication_level"
    dup_hdr_pos = grep(dup_hdr, hdr)
    if (length(dup_hdr_pos)) {
        len_rest = length(hdr)-dup_hdr_pos
        hdr = c(head(hdr, dup_hdr_pos-1L), paste0(dup_hdr, "_percent"), 
                    dup_hdr, tail(hdr, len_rest))
        idx1 = c(head(idx1, dup_hdr_pos), idx1[dup_hdr_pos]+1L, 
                    tail(idx1, len_rest))
        idx2 = c(head(idx2, dup_hdr_pos-1L), idx1[dup_hdr_pos], 
                    tail(idx2, len_rest+1L))
    }
    tables = lapply(seq_along(idx1), function(i) extract(idx1[i], idx2[i]))
    setattr(tables, 'names', hdr)
    tables = massage_table(tables)
    tables
}

massage_table <- function(x) {

    measure=NULL
    # widen basic statistics table and set column types properly
    y = x[["basic_statistics"]]
    y[, "measure" := gsub("[ ]+", "_", tolower(measure))]
    ans = dcast(y, group + sample_group + sample_name + pair ~ measure, 
            value.var = "value")
    cols = c("filtered_sequences", "percent_gc", "sequence_length", 
                "total_sequences")
    cols = intersect(cols, names(ans))
    ans[, (cols) := lapply(.SD, type.convert), .SDcols=cols]
    x[["basic_statistics"]] = ans
    x
}

# deprecated, but kept for now
split_sample <- function(sample_name) {
    splits = tstrsplit(sample_name, "_(?=[12]$)", perl=TRUE)
    if (length(splits) == 1L)
        splits = c(splits, list(1L))
    if (length(splits) != 2L)
        stop("'sample_name' argument should be of the form ", 
            "<samplename>_<pair> for paired end and <samplename> ", 
            "for single end Fastq files, with no '_' elsewhere.")
    splits[[2]] = as.integer(splits[[2]])
    splits
    # x[, c("sample_name", "pair") := splits]
    # setcolorder(x, c("sample_name", "pair", head(names(x), -2L)))
}
