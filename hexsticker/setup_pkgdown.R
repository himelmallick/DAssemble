# =============================================================================
# Reusable pkgdown setup for any R package
#
# Usage (from anywhere):
#   source("setup_pkgdown.R")
#   setup_pkgdown_site("~/Desktop/MyPkg", logo = "logo.png")      # local site
#   publish_pkgdown_site("~/Desktop/MyPkg")                        # GitHub Pages
#
# Requirements: install.packages(c("pkgdown", "usethis", "roxygen2", "yaml"))
# =============================================================================

library(usethis)
library(pkgdown)

# -----------------------------------------------------------------------------
# 1. Configure pkgdown, add the logo, and build the site locally
# -----------------------------------------------------------------------------
setup_pkgdown_site <- function(path      = ".",
                               logo      = NULL,          # path to a hex PNG
                               primary   = "#0E2A47",     # navbar / accent colour
                               bootswatch = "flatly",
                               font      = "Inter",
                               build     = TRUE,
                               preview   = interactive()) {
  path <- normalizePath(path, mustWork = TRUE)
  usethis::local_project(path, force = TRUE)
  pkg  <- desc::desc_get_field("Package", file = path)

  # creates _pkgdown.yml, ignores docs/ and adds Rbuildignore entries
  if (!file.exists(file.path(path, "_pkgdown.yml"))) usethis::use_pkgdown()

  # logo -> man/figures/logo.png (pkgdown picks it up automatically)
  if (!is.null(logo)) {
    usethis::use_logo(normalizePath(logo), retina = TRUE)

    # put the logo in the README title line (README.Rmd if present, else .md)
    readme <- file.path(path, c("README.Rmd", "README.md"))
    readme <- readme[file.exists(readme)][1]
    if (!is.na(readme)) {
      lines <- readLines(readme, warn = FALSE)
      h1 <- grep("^# ", lines)[1]
      if (!is.na(h1) && !grepl("man/figures/logo", lines[h1], fixed = TRUE)) {
        lines[h1] <- paste0(lines[h1], ' <img src="man/figures/logo.png" ',
                            'align="right" height="139" alt="" />')
        writeLines(lines, readme)
        cli::cli_alert_success("Added logo to {.file {basename(readme)}}")
      }
    }
  }

  # site URL from a GitHub remote, if there is one
  gh  <- tryCatch(usethis:::github_remote_list()$repo_spec[1],
                  error = function(e) NA_character_)
  url <- if (!is.na(gh)) {
    sprintf("https://%s.github.io/%s/",
            tolower(strsplit(gh, "/")[[1]][1]), strsplit(gh, "/")[[1]][2])
  }

  # reference index: every documented topic, the package overview last
  topics  <- sub("\\.Rd$", "", list.files(file.path(path, "man"), "\\.Rd$"))
  overview <- paste0(pkg, "-package")
  ref <- list(list(title = "Functions",
                   contents = as.list(setdiff(topics, overview))))
  if (overview %in% topics) {
    ref <- c(ref, list(list(title = "About the package",
                            contents = list(overview))))
  }

  cfg <- yaml::read_yaml(file.path(path, "_pkgdown.yml"))
  if (is.null(cfg)) cfg <- list()
  if (!is.null(url)) cfg$url <- url
  cfg$template <- list(
    bootstrap  = 5L,
    bootswatch = bootswatch,
    bslib = list(primary = primary,
                 base_font = list(google = font),
                 heading_font = list(google = font),
                 code_font = list(google = "JetBrains Mono"))
  )
  cfg$navbar$bg   <- "primary"   # dark navbar in the logo colour
  cfg$navbar$type <- "dark"
  # keep a hand-edited reference index if one already exists
  if (is.null(cfg$reference)) cfg$reference <- ref
  yaml::write_yaml(cfg, file.path(path, "_pkgdown.yml"))
  cli::cli_alert_success("Wrote {.file _pkgdown.yml}")

  if (build) {
    roxygen2::roxygenise(path)
    pkgdown::check_pkgdown(path)
    pkgdown::build_site(path, preview = preview)
  }
  invisible(path)
}

# -----------------------------------------------------------------------------
# 2. Publish to GitHub Pages (needs a GitHub remote + a PAT: usethis::gh_token_help())
#    Adds .github/workflows/pkgdown.yaml, creates the gh-pages branch, turns on
#    Pages, and writes the site URL into DESCRIPTION. Commit + push afterwards.
# -----------------------------------------------------------------------------
publish_pkgdown_site <- function(path = ".") {
  usethis::local_project(normalizePath(path), force = TRUE)
  usethis::use_pkgdown_github_pages()
}

# -----------------------------------------------------------------------------
# Example: DAssemble
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  setup_pkgdown_site("~/Desktop/DAssemble",
                     logo = "~/Desktop/DAssemble/hexsticker/logo_volcano.png")
}
