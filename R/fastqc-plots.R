#' @title Plot FastQC GC\% stats
#' @export
plot_gc_stats <- function(..., type=c("ggplot2", "rbokeh")) {

  ll = list(...)
  gc = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Basic Statistics")))
  if (is.null(names(ll)))
    setattr(gc, 'names', names(ll))
  gc = rbindlist(gc, idcol=TRUE)[`#Measure` == "%GC"]
  setnames(gc, "Value", "percent")[, percent := as.numeric(as.character(percent))]
  gc[, splits := findInterval(1:nrow(gc), seq(1, nrow(gc), by = 26L))]
  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(gc, aes(x=sample_name, y=percent)) + 
            geom_bar(stat="identity", aes(fill=group, group=sample_name)) + 
            ylim(0, 100) + 
            facet_wrap(splits ~ .id, scales="free_x", ncol=1L) + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("GC %") + 
             scale_fill_brewer(name = "Group", palette="Set1")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_total_sequence_stats <- function(..., type=c("ggplot2", "rbokeh")) {

  ll = list(...)
  ts = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Basic Statistics")))
  if (is.null(names(ll)))
    setattr(ts, 'names', names(ll))
  ts = rbindlist(ts, idcol=TRUE)[`#Measure` == "Total Sequences"]
  setnames(ts, "Value", "counts")[, counts  := as.numeric(as.character(counts))/1e6L]
  ts[, splits := findInterval(1:nrow(ts), seq(1, nrow(ts), by = 26L))]
  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(ts, aes(x=sample_name, y=counts)) + 
            geom_bar(stat="identity", aes(fill=group, group=sample_name)) + 
            facet_wrap(splits ~ .id, ncol=1L, scales="free_x") + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("Counts (in million)") + 
             scale_fill_brewer(name="Group", palette="Set1")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_dup_stats <- function(..., type=c("ggplot2", "rbokeh")) {

  ll = list(...)
  dup = lapply(ll, function(l) {
          rbindlist(lapply(l, function(y) { 
            v = y[["Sequence Duplication Levels"]][1]
            col = grep("[pP]ercentage", names(v))
            list(group=v$group, sample_group = v$sample_group, 
                 sample_name=v$sample_name, pair = v$pair, 
                 Duplication_Percent=as.numeric(names(v)[col+1L]))
          }))
        })
  if (is.null(names(ll)))
    setattr(dup, 'names', names(ll))
  dup = rbindlist(dup, idcol=TRUE)
  setnames(dup, "Duplication_Percent", "dup_percent")
  dup[, splits := findInterval(1:nrow(dup), seq(1, nrow(dup), by = 26L))]
  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(dup, aes(x=sample_name, y=dup_percent)) + 
            geom_bar(stat="identity", aes(fill=group, group=sample_name)) + 
            ylim(0, 100) + 
            facet_wrap(splits ~ .id, ncol=1L, scales="free_x") + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("Sample") + ylab("Duplication %") + 
             scale_fill_brewer(name="Group", palette="Set1")
    }, 
    rbokeh = {
      stop("Not yet implemented.")
    })
  pl
}

#' @title Plot FastQC per base sequence quality
#' @export
plot_sequence_quality <- function(..., type=c("ggplot2", "rbokeh")) {

  ll = list(...)
  seqn = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Per base sequence quality")))
  if (is.null(names(ll)))
    setattr(seqn, 'names', names(ll))
  seqn = rbindlist(seqn, idcol=TRUE)
  colnames = setdiff(names(seqn), c(".id", "group", "sample_group", "sample_name", "pair", "#Base"))
  seqn[, (colnames) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols=(colnames)]
  seqn[, pair := factor(pair)
    ][, `#Base` := factor(`#Base`, levels=unique(`#Base`))]
  val = rep(1:uniqueN(seqn[["sample_name"]]), each=nrow(seqn)/uniqueN(seqn[["sample_name"]]))
  seqn[, splits := findInterval(val, seq(1, nrow(seqn), by = 26L))]
  type = match.arg(type)
  switch(type, 
    ggplot2 = {
      if (!require(ggplot2))
        stop("Package 'ggplot2' is not available.")
      pl = ggplot(seqn, aes(x=`#Base`, colour=pair, group=sample_name)) + 
             geom_point(aes(y=Median)) + 
             geom_errorbar(aes(ymax=`90th Percentile`, ymin=`10th Percentile`)) + 
             geom_line(aes(y=Mean, group=sample_name)) + 
             facet_wrap(splits ~ .id, scales="free_x", ncol=1L) + 
             geom_hline(aes(yintercept=30L), colour="gray25", linetype="dotted") + 
             theme_bw() + 
             theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
             xlab("#Base") + ylab("Read quality") + 
             ylim(0L, max(seqn$Mean, na.rm=TRUE)+1L) + 
             scale_colour_brewer(name="Pair", palette="Set1")
    }, 
    rboken = {
      stop("Not yet implemented.")
    })
  pl
}
