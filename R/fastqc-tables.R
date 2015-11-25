#' Extract FastQC summary tables
#' @export
fastqc <- function(dir, groups=NULL) {
	# not exposing `id` argument here for now.
	files = list.files(dir, pattern="^fastqc_data\\.txt$", full.names=TRUE, recursive=TRUE)
	samples = gsub("_fastqc$", "", basename(dirname(files)))
	all_samples = lapply(seq_along(samples), function(i) fastqc_tables(samples[i], files[i], groups=groups))
	setattr(all_samples, 'class', 'fastqc')
	setattr(all_samples, 'names', samples)
}

fastqc_tables <- function(sample_name, fastqc_file, groups=NULL) {
	doc = readLines(fastqc_file)
	tables = extract_summary_tables(doc, sample_name, groups)
}

extract_summary_tables <- function(doc, sample_name, groups=NULL) {

	idx = grep("^>>", doc)

	idx1 = idx[c(TRUE, FALSE)]
	idx2 = idx[c(FALSE, TRUE)]

	hdr = gsub("^>>(.*)\t.*$", "\\1", doc[idx1])

	extract <- function(i, j) {
		if (j-i-2L <= 0L) {
			 message("Table at start index ", i, " seems to contain only one (header) row, therefore returning empty data.table.")
			 return(data.table:::null.data.table())
		}
		hdr_row = unlist(strsplit(doc[i+1L], "\t", fixed=TRUE))
		bdy_row = strsplit(doc[(i+2L):(j-1L)], "\t", fixed=TRUE)
		ans = setDT(transpose(bdy_row))
		setnames(ans, hdr_row)
		ans[, sample_name := sample_name
		  ][, c("sample_group", "pairs") := split_sample(sample_name)]
		if (is.data.table(groups)) {
			groups = copy(groups)[, sample_group := as.character(sample)]
			ans[groups, group := factor(i.group), on="sample_group"]
		} else {
			ans[, group := ""]
		}
		movecols = c("group", "sample_group", "sample_name", "pairs")
		setcolorder( ans, c( movecols, head(names(ans), -length(movecols)) ) )
	}
    if (is.data.table(groups) && !all(c("sample", "group") %in% names(groups)))
      stop("Columns 'sample' and 'group' should be present in argument 'groups'.")
	tables = lapply(seq_along(idx1), function(i) extract(idx1[i], idx2[i]))
	setattr(tables, 'names', hdr)
}

split_sample <- function(sample_name) {
	splits = tstrsplit(sample_name, "_(?=[12]$)", perl=TRUE)
	if (length(splits) == 1L)
		splits = c(splits, list(1L))
	if (length(splits) != 2L)
		stop("'sample_name' argument should be of the form <samplename>_<pair> for paired end and <samplename> for single end Fastq files, with no '_' elsewhere.")
	splits[[2]] = as.integer(splits[[2]])
	splits
	# x[, c("sample_name", "pairs") := splits]
	# setcolorder(x, c("sample_name", "pairs", head(names(x), -2L)))
}
