.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('v', as.character(packageVersion("ggfastqc")), 
        ', type vignette("ggfastqc-vignette", package="ggfastqc") to get started.')
    }
}
