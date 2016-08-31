.onLoad <- function(libname, pkgname) {
    ggfastqc_verbose = "FALSE"
    # global options
    opts = c(
              "ggfastqc_verbose" = ggfastqc_verbose
            )
    for (i in setdiff(names(opts), names(options())) ) {
        text = paste('options(', i, '=', opts[i], ')', sep="")
        eval(parse(text=text))
    }
    invisible()
}
