#' @title Plot FastQC GC\% stats
#' @export
plot_gc_stats <- function(..., type=c("ggplot2", "rbokeh"), group=NULL) {

  ll = list(...)
  gc = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Basic Statistics")))
  if (is.null(names(ll)))
    setattr(gc, 'names', names(ll))
  gc = rbindlist(gc, idcol=TRUE)[`#Measure` == "%GC"]
  gc[, c("sample_name", "pairs") := split_sample(sample)]
  setnames(gc, c("sample", "sample_name", "Value"), c("Sample (Pairs)", "Sample", "Percent"))
  gc[, Percent := as.numeric(as.character(Percent))]
  gc[, Sample  := as.character(Sample)]
  if (is.data.table(group)) {
    if (!all(c("Sample", "Condition") %in% names(group)))
      stop("Columns 'Sample' and 'Condition' should be present in data.table 'group'.")
    group = copy(group)[, Sample := as.character(Sample)][]
    gc[group, Condition := i.Condition, on="Sample"]
  } else gc[, Condition := ""]

  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(gc, aes(x=`Sample (Pairs)`, y=Percent))
      pl = pl + geom_bar(stat="identity", aes(fill=Condition, group=`Sample (Pairs)`))
      pl = pl + ylim(0, 100)
      pl = pl + facet_wrap(~ .id, ncol=1L, scales="free_x")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_total_sequence_stats <- function(..., type=c("ggplot2", "rbokeh"), group=NULL) {

  ll = list(...)
  ts = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Basic Statistics")))
  if (is.null(names(ll)))
    setattr(ts, 'names', names(ll))
  ts = rbindlist(ts, idcol=TRUE)[`#Measure` == "Total Sequences"]
  ts[, c("sample_name", "pairs") := split_sample(sample)]
  ts[, Value := as.numeric(as.character(Value))/1e6L]
  setnames(ts, c("sample", "sample_name", "Value"), c("Sample (Pairs)", "Sample", "Counts (*1e6)"))
  ts[, Sample := as.character(Sample)]
  if (is.data.table(group)) {
    if (!all(c("Sample", "Condition") %in% names(group)))
      stop("Columns 'Sample' and 'Condition' should be present in data.table 'group'.")
    group = copy(group)[, Sample := as.character(Sample)][]
    ts[group, Condition := i.Condition, on="Sample"]
  } else ts[, Condition := ""]

  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(ts, aes(x=`Sample (Pairs)`, y=`Counts (*1e6)`)) + 
            geom_bar(stat="identity", aes(fill=Condition, group=`Sample (Pairs)`)) + 
            facet_wrap(~ .id, ncol=1L, scales="free_x") + 
            theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("Counts (*1e6)") + scale_fill_brewer(palette="Set1")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_dup_stats <- function(..., type=c("ggplot2", "rbokeh"), group=NULL) {

  ll = list(...)
  dup = lapply(ll, function(l) {
          rbindlist(lapply(l, function(y) { 
            v = y[["Sequence Duplication Levels"]][1]
            list(sample=v$sample, Duplication_Percent=as.numeric(names(v)[3L]))
          }))
        })
  if (is.null(names(ll)))
    setattr(dup, 'names', names(ll))
  dup = rbindlist(dup, idcol=TRUE)[, c("sample_name", "pairs") := split_sample(sample)][]
  setnames(dup, c("sample", "sample_name", "Duplication_Percent"), c("Sample (Pairs)", "Sample", "Dup %"))
  dup[, Sample := as.character(Sample)]
  if (is.data.table(group)) {
    if (!all(c("Sample", "Condition") %in% names(group)))
      stop("Columns 'Sample' and 'Condition' should be present in data.table 'group'.")
    group = copy(group)[, Sample := as.character(Sample)][]
    dup[group, Condition := i.Condition, on="Sample"]
  } else dup[, Condition := ""]

  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(dup, aes(x=`Sample (Pairs)`, y=`Dup %`)) + 
            geom_bar(stat="identity", aes(fill=Condition, group=`Sample (Pairs)`)) + 
            ylim(0, 100) + 
            facet_wrap(~ .id, ncol=1L, scales="free_x") + 
             theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("Duplication %") + scale_fill_brewer(palette="Set1")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC per base sequence quality
#' @export
plot_sequence_quality <- function(..., type=c("ggplot2", "rbokeh"), group=NULL) {

  ll = list(...)
  seq = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Per base sequence quality")))
  if (is.null(names(ll)))
    setattr(seq, 'names', names(ll))
  seq = rbindlist(seq, idcol=TRUE)
  seq[, names(seq)[-(1:3)] := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols=names(seq)[-(1:3)]]
  seq[, c("sample_name", "pairs") := split_sample(sample)]
  setnames(seq, c("sample", "sample_name", "pairs"), c("Sample (Pairs)", "Sample", "Pairs"))
  seq[, Sample := as.character(Sample)
    ][, Pairs := factor(Pairs)
    ][, `#Base` := factor(`#Base`, levels=unique(`#Base`))]
  if (is.data.table(group)) {
    if (!all(c("Sample", "Condition") %in% names(group)))
      stop("Columns 'Sample' and 'Condition' should be present in data.table 'group'.")
    group = copy(group)[, Sample := as.character(Sample)][]
    seq[group, Condition := i.Condition, on="Sample"]
  } else seq[, Condition := ""]
  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(seq, aes(x=`#Base`, colour=Pairs, group=`Sample (Pairs)`)) + 
             geom_point(aes(y=Median)) + 
             geom_errorbar(aes(ymax=`90th Percentile`, ymin=`10th Percentile`)) + 
             geom_line(aes(y=Mean, group=`Sample (Pairs)`)) + 
             facet_wrap( ~ .id, scales="free_x", ncol=1L) + 
             geom_hline(aes(yintercept=30L), colour="gray25", linetype="dotted") + 
             theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("Read quality") + scale_fill_brewer(palette="Set1") + 
             ylim(0L, max(seq$Mean, na.rm=TRUE)+1L)
    }, 
    rboken = {
      stop("Not yet implemented.")
    })
  pl
}
