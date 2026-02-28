# Contributing to ensembleHTE

Thank you for considering contributing to ensembleHTE! This document provides guidelines for contributions.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear, descriptive title
- Steps to reproduce the bug
- Expected vs actual behavior
- Your R version and platform
- A minimal reproducible example (reprex)

### Suggesting Enhancements

Enhancement suggestions are welcome! Please:
- Check if the enhancement is already suggested
- Provide a clear use case
- Explain why it would be useful to users

### Pull Requests

1. Fork the repository
2. Create a new branch for your feature
3. Make your changes
4. Add tests for new functionality
5. Ensure `devtools::check()` passes
6. Update documentation
7. Submit a pull request

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "testthat", "usethis"))

# Clone and setup
git clone https://github.com/bfava/ensembleHTE.git
cd ensembleHTE
```

## Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use roxygen2 for documentation
- Write tests for all new functions
- Keep functions focused and modular

## Testing

All new code should include tests:

```r
# Run tests
devtools::test()

# Check test coverage
covr::package_coverage()
```

## Documentation

- Use roxygen2 for function documentation
- Include examples in `@examples`
- Update vignettes if adding major features
- Keep README.md up to date

## Code Review Process

All submissions require review. We aim to respond to pull requests within one week.

## Questions?

Feel free to open an issue for questions or discussion.
