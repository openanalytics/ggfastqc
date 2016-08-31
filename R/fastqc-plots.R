#' @title Plot FastQC GC\% stats
#' 
#' @description This function plots the \code{GC\%} of each of the samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{fastqc} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default), \code{"point"} 
#' and \code{"bar"}. \code{"jitter"} is only possible for 
#' \code{interactive = TRUE}, and is usually the preferred option since it 
#' provides the lowest ink ratio, and contains the least amount of clutter.
#' 
#' @examples
#' require(ggfastqc)
#' path <- system.file("tests/fastqc-sample", package="ggfastqc")
#' obj <- fastqc(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive = TRUE (plotly)
#' plot_gc_stats(sample = obj) # jitter
#' plot_gc_stats(sample = obj, geom = "point")
#' plot_gc_stats(sample = obj, geom = "bar")
#' 
#' # interactive = FALSE (ggplot2)
#' plot_gc_stats(sample = obj, interactive = FALSE) # jitter
#' plot_gc_stats(sample = obj, interactive = FALSE, geom = "point")
#' plot_gc_stats(sample = obj, interactive = FALSE, geom = "bar")
#' @seealso \code{\link{fastqc}} \code{\link{plot_dup_stats}}
#' \code{\link{plot_sequence_quality}} \code{\link{plot_total_sequence_stats}}
#' @export
plot_gc_stats <- function(..., interactive=TRUE, 
                    geom=c("jitter", "point", "bar")) {

    ll = list(...)
    gc = lapply(ll, function(l) {
            stopifnot(inherits(l, "fastqc"))
            ans = l[param == "basic_statistics", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(gc, 'names', paste0("fastqc_obj", seq_along(ll)))
    gc = rbindlist(gc, idcol=TRUE)
    cols = c("sample_name", "group", ".id")
    as_factor <- function(x) factor(x, levels=unique(x))
    gc[, (cols) := lapply(.SD, as_factor), .SDcols=cols
      ][, "splits" := findInterval(1:nrow(gc), seq(1, nrow(gc), by = 26L))]
    geom = match.arg(geom)

    aes = list(
            x = if(geom == "jitter") "group" else "sample_name", 
            y = "percent_gc", 
            sample_name = "sample_name", 
            fill = "group"
        )
    theme = list(
            xlab = if (geom == "jitter") "Groups" else "Sample", 
            ylab = "GC %", 
            title = "GC content among samples"
        )
    pl = switch(geom, 
            jitter = fastqc_jitter(gc, aes, theme, interactive), 
            point = fastqc_point(gc, aes, theme, interactive),
            bar = fastqc_bar(gc, aes, theme, interactive)
            # no need for stop(). match.arg() check above takes care of it.
        )
    pl
}

#' @title Plot FastQC total sequences stats
#' 
#' @description This function plots the total sequences across samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{fastqc} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default), \code{"point"} 
#' and \code{"bar"}. \code{"jitter"} is only possible for 
#' \code{interactive = TRUE}, and is usually the preferred option since it 
#' provides the lowest ink ratio, and contains the least amount of clutter.
#' 
#' @examples
#' require(ggfastqc)
#' path <- system.file("tests/fastqc-sample", package="ggfastqc")
#' obj <- fastqc(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive = TRUE (plotly)
#' plot_total_sequence_stats(sample = obj) # geom = "jitter" (default)
#' plot_total_sequence_stats(sample = obj, geom = "point")
#' plot_total_sequence_stats(sample = obj, geom = "bar")
#' 
#' # interactive = FALSE (ggplot2)
#' plot_total_sequence_stats(sample = obj, interactive = FALSE) # jitter
#' plot_total_sequence_stats(sample = obj, interactive = FALSE, geom = "point")
#' plot_total_sequence_stats(sample = obj, interactive = FALSE, geom = "bar")
#' @seealso \code{\link{fastqc}} \code{\link{plot_dup_stats}}
#' \code{\link{plot_sequence_quality}} \code{\link{plot_gc_stats}}
#' @export
plot_total_sequence_stats <- function(..., interactive=TRUE, 
                    geom=c("jitter", "point", "bar")) {

    ll = list(...)
    ts = lapply(ll, function(l) {
            stopifnot(inherits(l, "fastqc"))
            ans = l[param == "basic_statistics", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(ts, 'names', paste0("fastqc_obj", seq_along(ll)))
    ts = rbindlist(ts, idcol=TRUE)
    cols = c("sample_name", "group", ".id")
    as_factor <- function(x) factor(x, levels=unique(x))
    ts[, (cols) := lapply(.SD, as_factor), .SDcols=cols
      ][, "total_sequences"  := total_sequences/1e6L # in million reads
      ][, "splits" := findInterval(1:nrow(ts), seq(1, nrow(ts), by = 26L))]
    geom = match.arg(geom)

    aes = list(
            x = if(geom == "jitter") "group" else "sample_name", 
            y = "total_sequences", 
            sample_name = "sample_name", 
            fill = "group"
        )
    theme = list(
            xlab = if (geom == "jitter") "Groups" else "Sample", 
            ylab = "Sequence count (in million)", 
            title = "Sequence count among samples"
        )
    pl = switch(geom, 
            jitter = fastqc_jitter(ts, aes, theme, interactive), 
            point = fastqc_point(ts, aes, theme, interactive),
            bar = fastqc_bar(ts, aes, theme, interactive)
        )
    pl
}

#' @title Plot FastQC duplication percent stats
#' 
#' @description This function plots the duplication percent across samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{fastqc} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default), \code{"point"} 
#' and \code{"bar"}. \code{"jitter"} is only possible for 
#' \code{interactive = TRUE}, and is usually the preferred option since it 
#' provides the lowest ink ratio, and contains the least amount of clutter.
#' 
#' @examples
#' require(ggfastqc)
#' path <- system.file("tests/fastqc-sample", package="ggfastqc")
#' obj <- fastqc(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive = TRUE (plotly)
#' plot_dup_stats(sample = obj) # geom = "jitter" (default)
#' plot_dup_stats(sample = obj, geom = "point")
#' plot_dup_stats(sample = obj, geom = "bar")
#' 
#' # interactive = FALSE (ggplot2)
#' plot_dup_stats(sample = obj, interactive = FALSE) # jitter
#' plot_dup_stats(sample = obj, interactive = FALSE, geom = "point")
#' plot_dup_stats(sample = obj, interactive = FALSE, geom = "bar")
#' @seealso \code{\link{fastqc}} \code{\link{plot_total_sequence_stats}}
#' \code{\link{plot_sequence_quality}} \code{\link{plot_gc_stats}}
#' @export
plot_dup_stats <- function(..., interactive=TRUE, 
                    geom=c("jitter", "point", "bar")) {

    ll = list(...)
    dup = lapply(ll, function(l) {
            stopifnot(inherits(l, "fastqc"))
            ans = l[param == "sequence_duplication_level_percent", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(dup, 'names', paste0("fastqc_obj", seq_along(ll)))
    dup = rbindlist(dup, idcol=TRUE)
    cols = c("sample_name", "group", ".id")
    as_factor <- function(x) factor(x, levels=unique(x))
    dup[, (cols) := lapply(.SD, as_factor), .SDcols=cols
      ][, splits := findInterval(1:nrow(dup), seq(1, nrow(dup), by = 26L))]
    setnames(dup, "total_duplicate_percentage", "dup_percent")
    geom = match.arg(geom)

    aes = list(
            x = if(geom == "jitter") "group" else "sample_name", 
            y = "dup_percent", 
            sample_name = "sample_name", 
            fill = "group"
        )
    theme = list(
            xlab = if (geom == "jitter") "Groups" else "Sample", 
            ylab = "Sequence duplication (in %)", 
            title = "Sequence duplication among samples"
        )
    pl = switch(geom, 
            jitter = fastqc_jitter(dup, aes, theme, interactive), 
            point = fastqc_point(dup, aes, theme, interactive),
            bar = fastqc_bar(dup, aes, theme, interactive)
        )
    pl
}

#' @title Plot FastQC per base sequence quality
#' 
#' @description This function plots the per base sequence quality across 
#' samples. If info regarding the \code{group} each \code{sample} belongs to 
#' is also available, then the generated plot will take that into account to 
#' colour / facet accordingly.
#' 
#' @param \dots The set of \code{fastqc} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Only possible value is \code{"line"}.
#' 
#' @examples
#' require(ggfastqc)
#' path <- system.file("tests/fastqc-sample", package="ggfastqc")
#' obj <- fastqc(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive = TRUE (plotly)
#' plot_sequence_quality_stats(sample = obj)
#' 
#' # interactive = FALSE (ggplot2)
#' plot_sequence_quality_stats(sample = obj, interactive = FALSE)
#' @seealso \code{\link{fastqc}} \code{\link{plot_dup_stats}}
#' \code{\link{plot_total_sequence_stats}} \code{\link{plot_gc_stats}}
#' @export
plot_sequence_quality <- function(..., interactive=TRUE, geom=c("line")) {

    ll = list(...)
    seqn = lapply(ll, function(l) {
            stopifnot(inherits(l, "fastqc"))
            ans = l[param == "per_base_sequence_quality", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(seqn, 'names', paste0("fastqc_obj", seq_along(ll)))
    seqn = rbindlist(seqn, idcol=TRUE)
    seqn[, base_int := as.integer(gsub("_.*$", "", base))]
    cols = c("sample_name", "group", ".id", "pair")
    as_factor <- function(x) factor(x, levels=unique(x))
    seqn[, (cols) := lapply(.SD, as_factor), .SDcols=cols]
    val = rep(1:uniqueN(seqn[["sample_name"]]), 
            each=nrow(seqn)/uniqueN(seqn[["sample_name"]]))
    seqn[, splits := findInterval(val, seq(1, nrow(seqn), by = 26L))]
    geom = match.arg(geom)

    aes = list(
            x = "base_int",  
            y = "mean", 
            colour = "pair", 
            group = "sample_name", 
            ymax = "`90th_percentile`", 
            ymin = "`10th_percentile`", 
            sample_name = "sample_name", 
            base = "base"
        )
    theme = list(
            xlab = "# Base", 
            ylab = "Read quality", 
            title = "Per base sequence quality"
        )
    pl = switch(geom, 
            line = fastqc_errorbar(seqn, aes, theme, interactive), 
        )
    pl
}

## internal geoms ----------------------- 

fastqc_jitter <- function(dt, aes, theme, interactive=TRUE) {

    if (!interactive) {
        warn = paste("When geom='jitter', only interactive", 
                "plots are possible. Setting interactive=TRUE")
        warning(warn)
        interactive=TRUE
    }
    if (!requireNamespace("plotly"))
        stop("Package 'plotly' is not available.")
    p = ggplot(dt, aes_string(x=aes$x, y=aes$y, sample=aes$sample_name)) + 
        geom_point(position = position_jitter(), aes_string(fill=aes$fill), 
            shape=21L, size=4L) + 
        ylim(0L, 100L) + xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        facet_wrap(~ .id, scales="free_y", ncol=1L) + 
        theme(axis.text.x = element_text(angle = 45L, hjust = 1L), 
            legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="none") + 
        ggtitle(theme$title)

    plotly::ggplotly(p, tooltip = c("sample", "y"))
}

fastqc_point <- function(dt, aes, theme, interactive=TRUE) {
    p = ggplot(dt, aes_string(x=aes$x, y=aes$y)) + 
        geom_point(aes_string(fill=aes$fill), shape=21L, size=4L) + 
        ylim(0L, 100L) + xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        (if (max(dt$splits) == 1L)
            facet_wrap(~ .id, scales="free_y", ncol=1L)
        else  facet_wrap(splits ~ .id, scales="free_x", ncol=1L)) + 
        theme(axis.text.x = element_text(angle = 45L, hjust = 1L), 
            legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p, tooltip = c("x", "y"))
    }
    p
}

fastqc_bar <- function(dt, aes, theme, interactive=TRUE) {
    p = ggplot(dt, aes_string(x=aes$x, y=aes$y)) + 
        geom_bar(stat="identity", aes_string(fill=aes$fill)) + 
        ylim(0L, 100L) + xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        (if (max(dt$splits) == 1L)
            facet_wrap(~ .id, scales="free_y", ncol=1L)
        else  facet_wrap(splits ~ .id, scales="free_x", ncol=1L)) + 
        theme(axis.text.x = element_text(angle = 45L, hjust = 1L), 
            legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p, tooltip = c("x", "y"))
    }
    p
}

fastqc_errorbar <- function(dt, aes, theme, interactive=TRUE) {
    p = ggplot(dt, aes_string(x=aes$x, colour=aes$colour, group=aes$group, 
            sample_name=aes$sample_name, percentile_10=aes$ymin, 
            percentile_90=aes$ymax, base=aes$base)) + 
        geom_point(aes_string(y=aes$y)) + 
        geom_errorbar(aes_string(ymax=aes$ymax, ymin=aes$ymin))+
        geom_line(aes_string(y=aes$y, group=aes$group)) + 
        (if (max(dt$splits) == 1L)
            facet_wrap(~ .id, scales="free_y", ncol=1L)
        else  facet_wrap(splits ~ .id, scales="free_x", ncol=1L)) + 
        geom_hline(aes(yintercept=30L), colour="gray25", linetype="dotted") + 
        ylim(0, max(dt$mean, na.rm=TRUE)+1L) + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 45L, hjust = 1L), 
            legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06)) + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        # no error bars in plotly yet.. tooltip includes base#, 
        # 10th and 90th percentile
        p = plotly::ggplotly(p, tooltip = 
                c("sample_name", "base", "percentile_10", "percentile_90"))
    }
    p
}
