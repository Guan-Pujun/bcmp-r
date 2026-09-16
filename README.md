<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/bcmp-github-dark.svg">
    <img src="docs/assets/bcmp-github-light.svg" alt="BCMP — Batch-Constraint Manifold Partition" width="460">
  </picture>
</p>

<br>

<p align="center">
  <a href="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/ci.yml"><img src="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/ci.yml/badge.svg?branch=main" alt="Tests"></a>
  <a href="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/package.yml"><img src="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/package.yml/badge.svg?branch=main" alt="Build"></a>
  <a href="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/docs.yml"><img src="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/docs.yml/badge.svg?branch=main" alt="Documentation"></a>
  <br>
  <a href="https://codecov.io/github/Guan-Pujun/bcmp-r"><img src="https://codecov.io/github/Guan-Pujun/bcmp-r/graph/badge.svg?token=YO04KDP79U" alt="Codecov"></a>
  <a href="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/R-CMD-check.yml"><img src="https://github.com/Guan-Pujun/bcmp-r/actions/workflows/R-CMD-check.yml/badge.svg?branch=main" alt="R CMD check"></a>
  <a href="https://www.bestpractices.dev/projects/14679"><img src="https://www.bestpractices.dev/projects/14679/badge" alt="OpenSSF Best Practices"></a>
</p>

<br>

# BCMP<sup>R</sup>

Batch-Constraint Manifold Partition (BCMP) derives domain labels from the
unintegrated geometry of single-cell data. These labels can be used as proxy
population labels when evaluating batch-effect removal by candidate integration
methods.

## Installation

```r
install.packages("bcmp")
```

See [Installation](https://guan-pujun.github.io/bcmp-r/installation.html) for details.

## Quick start

For a Seurat object containing raw counts and a `batch` metadata column:

```r
library(bcmp)

result <- bcmp(seurat, batch_key = "batch")
result$object[["bcmp_domain"]]
```

`bcmp()` also supports SingleCellExperiment objects.

## Manual

The [BCMP<sup>R</sup> manual](https://guan-pujun.github.io/bcmp-r/)
describes installation, the public API, returned results, and the vignette.

## Feedback and contributions

Report bugs or request enhancements through the [GitHub issue tracker](https://github.com/Guan-Pujun/bcmp-r/issues).
For proposed changes, see [CONTRIBUTING.md](CONTRIBUTING.md).
Please report suspected security vulnerabilities through the private process in
[SECURITY.md](SECURITY.md), not in a public issue.

## Citation

Citation information will be added when available.

## License

BCMP is distributed under the [GPL-3.0-or-later license](LICENSE).
