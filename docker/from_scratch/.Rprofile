## Set CRAN Mirror:
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org/"
    options(repos = r)
})


if (interactive()) {
    .Last <- function() try(savehistory("~/.Rhistory"))
}
options(menu.graphics = FALSE)

.libPaths(c("/usr/lib/R/library/", "/usr/lib/R/site-library/"))