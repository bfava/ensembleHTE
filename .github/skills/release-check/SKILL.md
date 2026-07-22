---
name: release-check
description: "Pre-push release checklist for ensembleHTE. Use when: updating GitHub, pushing to main, preparing a release, before git push, before a new version, sanity-check before publishing. Runs document/test/check, verifies vignette fits, NEWS, version bump, and clean git status."
---

# Release Check (pre-GitHub-push)

Run through this checklist any time the user wants to push the package to GitHub or prepare a release. Stop at the first failure and report it — do not proceed silently.

## When to Use

- User says: "update GitHub", "push to main", "prepare release", "ready to publish", "bump version".
- Before tagging a new version or submitting to CRAN.

## Procedure

Execute each step in order. After each step, report PASS / FAIL with a one-line summary.

### 1. Clean working tree check
```bash
git status --short
git log origin/main..HEAD --oneline
```
Flag any untracked files outside `data-raw/`, `vignettes/precomputed/`, or `.Rcheck/`. Note unpushed commits.

### 2. Regenerate docs
```r
devtools::document()
```
If `man/*.Rd` or `NAMESPACE` changed, those changes must be staged. Hand-edits to either are a hard fail.

### 3. Style / lint (best-effort, non-blocking)
```r
if (requireNamespace("lintr", quietly = TRUE)) lintr::lint_package()
```

### 4. Run full test suite (not just CRAN-safe subset)
```bash
NOT_CRAN=true Rscript -e "devtools::test()"
```
All tests must pass. `skip_on_cran()` tests SHOULD run locally — investigate any unexpected skips.

### 5. R CMD check — must be 0 / 0 / 0
```r
devtools::check()
```
Required: 0 errors, 0 warnings, 0 notes. Update [cran-comments.md](../../../cran-comments.md) if the count changes.

### 6. Vignettes are pre-computed and current
- Confirm [vignettes/precomputed/](../../../vignettes/precomputed/) `.rds` files exist.
- If `R/ensemble_hte.R` or `R/ensemble_pred.R` public API changed since the last vignette refresh, OR if [data/microcredit.rda](../../../data/microcredit.rda) was rebuilt, re-run:
  ```r
  source("data-raw/precompute_vignette_fits.R")
  source("data-raw/precompute_walkthrough.R")
  ```
- Build vignettes to confirm they render:
  ```r
  devtools::build_vignettes()
  ```

### 7. Version & NEWS
- Inspect [DESCRIPTION](../../../DESCRIPTION) `Version:` field.
- If this push is a release (not a bugfix-in-progress), ask the user whether to bump (`usethis::use_version("patch" | "minor" | "major")`).
- Confirm [NEWS.md](../../../NEWS.md) has an entry under a heading matching the current version. If missing, prompt the user for bullet points and add them.
- If version bumped, update [inst/update_status.dcf](../../../inst/update_status.dcf) so `ensemble_news()` reflects the release.

### 8. pkgdown sanity (if site config touched)
Only if `_pkgdown.yml`, `README.md`, or any `vignettes/articles/**` changed:
```r
pkgdown::build_site(preview = FALSE)
```
The `pkgdown.yaml` GitHub Action will rebuild on push, so this is a local preview only.

### 9. Final git status
```bash
git status --short
```
Everything regenerated in steps 2 and 6 must be staged (or intentionally excluded). The `ensembleHTE.Rcheck/` directory must NOT be staged — it is local check output only.

### 10. Summary report
Print a markdown table with one row per step (PASS / FAIL / SKIPPED + note). Do NOT run `git push` automatically — surface the table and let the user push.

### 11. Publish the release on GitHub (tag is NOT enough)

A pushed git tag does **not** create a GitHub Release, and the repo's landing page / "Releases" sidebar shows the latest **Release**, not the latest tag. If you stop at `git push origin vX.Y.Z`, GitHub will keep displaying the previous version.

After the tag is pushed, always create the matching Release:

```bash
# Verify what GitHub currently shows as "Latest"
gh release list --limit 5

# Create the Release for the new tag, with notes from the NEWS.md section
gh release create vX.Y.Z --title "ensembleHTE X.Y.Z" --notes "<paste the NEWS.md X.Y.Z bullets>"

# Confirm the new version is now marked "Latest"
gh release list --limit 5
```

The release is only complete when `gh release list` shows the new version tagged `Latest`. `gh` is authenticated on this machine (`gh auth status`); if it is not, tell the user to run `gh auth login` themselves — never do it silently.

## Pitfalls

- Pushing a git tag is NOT a release. The GitHub landing page shows the latest **GitHub Release**, so you must run `gh release create vX.Y.Z` (step 11) or the site keeps showing the old version. Verify with `gh release list` that the new version is `Latest`.
- Running `devtools::test()` without `NOT_CRAN=true` skips the heavy tests that exercise `grf` and ensemble paths.
- Forgetting to re-run `precompute_vignette_fits.R` after an API change → vignettes break on the next pkgdown build.
- Staging `man/*.Rd` edits without re-running `devtools::document()` produces a doc/code mismatch.
- `_pkgdown.yml` builds an online site — do not commit secrets or local file paths into examples.
