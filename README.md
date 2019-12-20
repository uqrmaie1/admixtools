
<!-- README.md is generated from README.Rmd. Please edit that file -->

# admixtools

Lightweight R implementations of some
[AdmixTools](https://github.com/DReichLab/AdmixTools) programs.

## Installation

To install `admixtools`, first you need `devtools`:

``` r
install.packages("devtools")
```

Then install `admixtools`, located here on
O2:

``` r
pkg = 'file:///n/groups/reich/robert/projects/admixprograms/admixtools_0.1.0.tar.gz'
devtools::install_url(pkg, dependencies='Imports', upgrade='never')
```

Load it like this:

``` r
library("admixtools")
```

## Usage

If you already have precomputed f2-statistics, you can run estimate
admixture weights like this.

``` r
target = 'Denisova.DG'
left = c('Altai_Neanderthal.DG', 'Vindija.DG')
right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
```

``` r
f2_dir = '/n/groups/reich/robert/projects/admixprograms/f2blocks_v41.1/'
qpadm(target, left, right, f2_dir = f2_dir)
```

Or you can use these f2-statistics to fit an admixturegraph to the data.

``` r
qpg_results = qpgraph(graph1, f2_dir = f2_dir)
```

More examples are in the [`Get started`](articles/admixtools.html)
section.
