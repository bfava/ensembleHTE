# Check for ensembleHTE Updates

Checks whether a newer version of ensembleHTE is available on GitHub and
reports update urgency. Can also display the full changelog.

## Usage

``` r
ensemble_news(show_changelog = FALSE)
```

## Arguments

- show_changelog:

  Logical. If `TRUE`, displays the full `NEWS.md` changelog instead of
  checking for updates. Default is `FALSE`.

## Value

Invisibly returns `NULL`. Called for its side effect of printing update
status or changelog to the console.

## Details

The function fetches a small status file from the public GitHub
repository to determine whether an update is needed. It scans *all*
versions between the user's installed version and the latest release. If
*any* of those missed versions is marked `critical`, the user is told a
critical update is needed (even if the latest version itself is only
`recommended`).

Update statuses (per version):

- **critical**: A bug fix or breaking change. Users should update
  immediately.

- **recommended**: New functionality or improvements worth updating for.

- **current**: No action needed (used for the initial release).

If the network is unavailable or the repository cannot be reached, the
function reports the installed version and suggests viewing the local
changelog.

## For the package maintainer

To signal updates to users, edit `inst/update_status.dcf` in the
repository before pushing a new release. The file uses R's DCF format
with one record per version, separated by blank lines. Newest version
first:

    Version: 0.3.0
    Status: recommended
    Message: New feature: CLAN improvements.

    Version: 0.2.1
    Status: critical
    Message: Critical bug fix in gates() subsetting.

    Version: 0.2.0
    Status: recommended
    Message: Added gavs_restricted().

    Version: 0.1.0
    Status: current
    Message: Initial release.

A user on v0.2.0 running `ensemble_news()` would see "CRITICAL UPDATE
NEEDED" because v0.2.1 (which they missed) is critical, even though the
latest v0.3.0 is only recommended.

## Examples

``` r
if (FALSE) { # \dontrun{
# Check for updates
ensemble_news()

# View full changelog
ensemble_news(show_changelog = TRUE)
} # }
```
