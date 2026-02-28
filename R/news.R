#' Check for ensembleHTE Updates
#'
#' @description
#' Checks whether a newer version of ensembleHTE is available on GitHub and
#' reports update urgency. Can also display the full changelog.
#'
#' @param show_changelog Logical. If \code{TRUE}, displays the full \code{NEWS.md}
#'   changelog instead of checking for updates. Default is \code{FALSE}.
#'
#' @details
#' The function fetches a small status file from the public GitHub repository to
#' determine whether an update is needed. It scans \emph{all} versions between
#' the user's installed version and the latest release. If \emph{any} of those
#' missed versions is marked \code{critical}, the user is told a critical update
#' is needed (even if the latest version itself is only \code{recommended}).
#'
#' Update statuses (per version):
#' \itemize{
#'   \item \strong{critical}: A bug fix or breaking change. Users should update immediately.
#'   \item \strong{recommended}: New functionality or improvements worth updating for.
#'   \item \strong{current}: No action needed (used for the initial release).
#' }
#'
#' If the network is unavailable or the repository cannot be reached, the function
#' reports the installed version and suggests viewing the local changelog.
#'
#' @section For the package maintainer:
#' To signal updates to users, edit \code{inst/update_status.dcf} in the repository
#' before pushing a new release. The file uses R's DCF format with one record per
#' version, separated by blank lines. Newest version first:
#' \preformatted{
#' Version: 0.3.0
#' Status: recommended
#' Message: New feature: CLAN improvements.
#'
#' Version: 0.2.1
#' Status: critical
#' Message: Critical bug fix in gates() subsetting.
#'
#' Version: 0.2.0
#' Status: recommended
#' Message: Added gavs_restricted().
#'
#' Version: 0.1.0
#' Status: current
#' Message: Initial release.
#' }
#'
#' A user on v0.2.0 running \code{ensemble_news()} would see
#' "CRITICAL UPDATE NEEDED" because v0.2.1 (which they missed) is critical,
#' even though the latest v0.3.0 is only recommended.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of printing
#'   update status or changelog to the console.
#'
#' @examples
#' \dontrun{
#' # Check for updates
#' ensemble_news()
#'
#' # View full changelog
#' ensemble_news(show_changelog = TRUE)
#' }
#'
#' @export
ensemble_news <- function(show_changelog = FALSE) {

  installed_ver <- as.character(utils::packageVersion("ensembleHTE"))

  # --- Show changelog mode ---
  if (show_changelog) {
    news_file <- system.file("NEWS.md", package = "ensembleHTE")
    if (nzchar(news_file) && file.exists(news_file)) {
      cat(readLines(news_file, warn = FALSE), sep = "\n")
    } else {
      cat("No changelog found. See https://github.com/bfava/ensembleHTE for updates.\n")
    }
    return(invisible(NULL))
  }

  # --- Update check mode ---
  update_info <- .fetch_update_status()

  if (is.null(update_info) || !"Version" %in% colnames(update_info)) {
    cat("ensembleHTE v", installed_ver, "\n", sep = "")
    cat("Unable to check for updates (no network or repository unreachable).\n")
    cat('Run ensemble_news(show_changelog = TRUE) to view the local changelog.\n')
    return(invisible(NULL))
  }

  # Parse all version records
  versions  <- trimws(update_info[, "Version"])
  statuses  <- tolower(trimws(update_info[, "Status"]))
  messages  <- if ("Message" %in% colnames(update_info)) trimws(update_info[, "Message"]) else rep("", length(versions))
  latest_ver <- versions[1]

  # Find versions the user has missed (strictly newer than installed)
  missed <- vapply(versions, function(v) utils::compareVersion(v, installed_ver) > 0, logical(1))

  if (!any(missed)) {
    cat("ensembleHTE v", installed_ver, " -- up to date.\n", sep = "")
    return(invisible(NULL))
  }

  # Check if any missed version is critical
  has_critical <- any(statuses[missed] == "critical")

  if (has_critical) {
    # Find the critical version(s) for the message
    crit_idx <- which(missed & statuses == "critical")
    cat("!! CRITICAL UPDATE NEEDED !!\n")
    cat("   Installed: v", installed_ver, "\n", sep = "")
    cat("   Latest:    v", latest_ver, "\n", sep = "")
    cat("\n   Critical fix(es) you are missing:\n")
    for (i in crit_idx) {
      cat("   - v", versions[i], ": ", messages[i], "\n", sep = "")
    }
  } else {
    cat("Update recommended.\n")
    cat("   Installed: v", installed_ver, "\n", sep = "")
    cat("   Latest:    v", latest_ver, "\n", sep = "")
    # Show the latest version's message
    if (nzchar(messages[1])) {
      cat("   ", messages[1], "\n", sep = "")
    }
  }

  cat("\n   Update with:\n")
  cat('   devtools::install_github("bfava/ensembleHTE")\n')

  invisible(NULL)
}


#' Fetch update status from GitHub
#'
#' @description
#' Internal helper that fetches \code{inst/update_status.dcf} from the public
#' GitHub repository.
#'
#' @return A matrix from \code{read.dcf()}, or \code{NULL} on failure.
#' @keywords internal
.fetch_update_status <- function() {
  tryCatch({
    url_str <- "https://raw.githubusercontent.com/bfava/ensembleHTE/main/inst/update_status.dcf"
    con <- url(url_str)
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    if (length(lines) == 0) return(NULL)
    read.dcf(textConnection(paste(lines, collapse = "\n")))
  }, error = function(e) NULL)
}
