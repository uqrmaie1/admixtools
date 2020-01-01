
<!-- README.md is generated from README.Rmd. Please edit that file -->

# admixtools

Lightweight R implementations of some
[AdmixTools](https://github.com/DReichLab/AdmixTools) programs.

## Installation

To install and load `admixtools`, follow these steps:

``` r
install.packages(c('tidyverse', 'devtools'))
pkg = 'file:///n/groups/reich/robert/projects/admixprograms/admixtools_0.1.0.tar.gz'
devtools::install_url(pkg, dependencies='Imports', upgrade='never')
library('admixtools')
```

## Usage

If you already have precomputed f2-statistics, you can fit an admixture
graph like this:

``` r
f2_dir = '/n/groups/reich/robert/projects/admixprograms/f2blocks_v42.1/'
fit = qpgraph(graph1, f2_dir = f2_dir)
plot_graph(fit$edges)
```

![example graph](man/figures/graph1.png)

<br>

You can also use these f2-statistics to estimate admixture weights:

``` r
target = 'Denisova.DG'
left = c('Altai_Neanderthal.DG', 'Vindija.DG')
right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
qpadm(target, left, right, f2_dir = f2_dir)
```

More documentation
[here](https://uqrmaie1.github.io/admixtools/articles/admixtools.html).
