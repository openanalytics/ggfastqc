#' @title Plot FastQC GC\% stats
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
    pl = switch(geom, 
        jitter = {
            if (!interactive) {
                warn = paste("When geom='jitter', only interactive", 
                        "plots are possible. Setting interactive=TRUE")
                warning(warn)
                interactive=TRUE
            }
            if (!requireNamespace("plotly"))
                stop("Package 'plotly' is not available.")
            p = ggplot(gc, aes(x=group, y=percent_gc, sample=sample_name)) + 
                geom_point(position = position_jitter(), aes(fill=group), 
                    shape=21L, size=4L) + 
                ylim(0, 100) + xlab("Sample") + ylab("GC %") + theme_bw() + 
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
                ggtitle("GC content distribution among samples")
            plotly::ggplotly(p, tooltip = c("sample", "y"))
        }, 
        point = {
            p = ggplot(gc, aes(x=sample_name, y=percent_gc)) + 
                geom_point(aes(fill=group), shape=21L, size=4L) + 
                ylim(0, 100) + xlab("Sample") + ylab("GC %") + theme_bw() + 
                (if (max(gc$splits) == 1L)
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
                ggtitle("GC content distribution among samples")
            if (interactive)
                p = plotly::ggplotly(p, tooltip = c("x", "y"))
            p
        }, 
        bar = {
            p = ggplot(gc, aes(x=sample_name, y=percent_gc)) + 
                geom_bar(stat="identity", aes(fill=group)) + 
                ylim(0, 100) + xlab("Sample") + ylab("GC %") + theme_bw() + 
                (if (max(gc$splits) == 1L)
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
                ggtitle("GC content distribution among samples")
            if (interactive)
                p = plotly::ggplotly(p, tooltip = c("x", "y"))
            p
        }, 
    )
    pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_total_sequence_stats <- function(..., type=c("ggplot2", "plotly")) {

    ll = list(...)
    ts = lapply(ll, function(l) rbindlist(lapply(l, `[[`, "Basic Statistics")))
    if (is.null(names(ll)))
        setattr(ts, 'names', names(ll))
    ts = rbindlist(ts, idcol=TRUE)[`#Measure` == "Total Sequences"]
    setnames(ts, "Value", "counts")
    ts[, "counts"  := as.numeric(as.character(counts))/1e6L
      ][, "splits" := findInterval(1:nrow(ts), seq(1, nrow(ts), by = 26L))]
    type = match.arg(type)
    if (!requireNamespace(ggplot2))
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
    if (type == "plotly") {
        if (!requireNamespace(plotly))
            stop("Package 'plotly' is not available.")
        pl = plotly::ggplotly(pl)
    }
    pl
}

#' @title Plot FastQC total sequences stats
#' @export
plot_dup_stats <- function(..., type=c("ggplot2", "plotly")) {

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
    if (!requireNamespace(ggplot2))
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
    if (type == "plotly") {
        if (!requireNamespace(plotly))
            stop("Package 'plotly' is not available.")
        pl = plotly::ggplotly(pl)
    }
    pl
}

#' @title Plot FastQC per base sequence quality
#' @export
plot_sequence_quality <- function(..., type=c("ggplot2", "plotly")) {

    ll = list(...)
    seqn = lapply(ll, function(l) 
            rbindlist(lapply(l, `[[`, "Per base sequence quality")))
    if (is.null(names(ll)))
        setattr(seqn, 'names', names(ll))
    seqn = rbindlist(seqn, idcol=TRUE)
    cols = c(".id", "group", "sample_group", "sample_name", "pair", "#Base")
    cols = setdiff(names(seqn), cols)
    fac2num <- function(x) as.numeric(as.character(x))
    seqn[, (cols) := lapply(.SD, fac2num, .SDcols=cols]
    seqn[, pair := factor(pair)
        ][, `#Base` := factor(`#Base`, levels=unique(`#Base`))]
    val = rep(1:uniqueN(seqn[["sample_name"]]), 
            each=nrow(seqn)/uniqueN(seqn[["sample_name"]]))
    seqn[, splits := findInterval(val, seq(1, nrow(seqn), by = 26L))]
    type = match.arg(type)
    if (!requireNamespace(ggplot2))
        stop("Package 'ggplot2' is not available.")
    pl = ggplot(seqn, aes(x=`#Base`, colour=pair, group=sample_name)) + 
            geom_point(aes(y=Median)) + 
            geom_errorbar(aes(ymax=`90th Percentile`, ymin=`10th Percentile`))+
            geom_line(aes(y=Mean, group=sample_name)) + 
            facet_wrap(splits ~ .id, scales="free_x", ncol=1L) + 
            geom_hline(aes(yintercept=30L), colour="gray25", linetype="dotted")+
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45L, hjust = 1, size=14), 
                axis.text.y = element_text(size = 14),  
                text = element_text(size=18), 
                panel.grid.major = element_blank()) + 
            xlab("#Base") + ylab("Read quality") + 
            ylim(0L, max(seqn$Mean, na.rm=TRUE)+1L) + 
            scale_colour_brewer(name="Pair", palette="Set1")
    if (type == "plotly") {
        if (!requireNamespace(plotly))
            stop("Package 'plotly' is not available.")
        pl = plotly::ggplotly(pl)
    }
    pl
}
