
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

Or, if you don’t want to install it, the following might work instead if
you use R 3.5 on
O2:

``` r
.libPaths(c(.libPaths(), '/home/rm360/R/x86_64-pc-linux-gnu-library/3.5/'))
library('admixtools')
```

## Usage

If you already have precomputed f2-statistics, you can fit an admixture
graph like this:

``` r
f2_dir = '/n/groups/reich/robert/projects/admixprograms/f2blocks_v42.1/'
fit = qpgraph(f2_dir, example_graph)
plot_graph(fit$edges)
```

![example graph](man/figures/graph1.png)

<br>

You can also use the f2-statistics to estimate admixture weights:

``` r
target = 'Denisova.DG'
left = c('Altai_Neanderthal.DG', 'Vindija.DG')
right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
```

``` r
qpadm(f2_dir, target, left, right)
```

    #> # A tibble: 2 x 4
    #>   target      left                 weight    se
    #>   <chr>       <chr>                 <dbl> <dbl>
    #> 1 Denisova.DG Altai_Neanderthal.DG   43.4  17.7
    #> 2 Denisova.DG Vindija.DG            -42.4  17.7

More documentation
[here](https://uqrmaie1.github.io/admixtools/articles/admixtools.html).
