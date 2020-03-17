
<!-- README.md is generated from README.Rmd. Please edit that file --->

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

    #> ℹ Computing f4 stats...
    #> ℹ Computing admixture weights...
    #> ℹ Computing standard errors...
    #> ℹ Computing number of admixture waves...
    #> $weights
    #> # A tibble: 2 x 5
    #>   target      left                 weight    se     z
    #>   <chr>       <chr>                 <dbl> <dbl> <dbl>
    #> 1 Denisova.DG Altai_Neanderthal.DG   50.2  24.0  2.10
    #> 2 Denisova.DG Vindija.DG            -49.2  24.0 -2.06
    #> 
    #> $f4
    #> # A tibble: 36 x 9
    #>    pop1    pop2       pop3    pop4            est      se     z         p weight
    #>    <chr>   <chr>      <chr>   <chr>         <dbl>   <dbl> <dbl>     <dbl>  <dbl>
    #>  1 Deniso… Altai_Nea… Chimp.… Mbuti.DG    0.0101  2.82e-4  35.6 4.79e-278   50.2
    #>  2 Deniso… fit        Chimp.… Mbuti.DG   -0.00616 3.91e-4 -15.7 9.78e- 56   NA  
    #>  3 Deniso… Vindija.DG Chimp.… Mbuti.DG    0.0102  2.90e-4  35.1 1.28e-269  -49.2
    #>  4 Deniso… Altai_Nea… Chimp.… Russia_Us…  0.0118  3.48e-4  34.0 4.75e-253   50.2
    #>  5 Deniso… fit        Chimp.… Russia_Us…  0.0144  8.93e-4  16.1 1.64e- 58   NA  
    #>  6 Deniso… Vindija.DG Chimp.… Russia_Us…  0.0122  3.53e-4  34.5 4.40e-260  -49.2
    #>  7 Deniso… Altai_Nea… Chimp.… Switzerla…  0.0116  3.65e-4  31.9 1.88e-223   50.2
    #>  8 Deniso… fit        Chimp.… Switzerla…  0.0846  4.96e-4 171.  0.          NA  
    #>  9 Deniso… Vindija.DG Chimp.… Switzerla…  0.0120  3.76e-4  31.9 1.14e-222  -49.2
    #> 10 Deniso… Altai_Nea… Mbuti.… Chimp.REF  -0.0101  2.82e-4 -35.6 4.79e-278   50.2
    #> # … with 26 more rows
    #> 
    #> $rankdrop
    #> # A tibble: 1 x 7
    #>   f4rank   dof chisq      p dofdiff chisqdiff p_nested
    #> *  <int> <int> <dbl>  <dbl>   <int>     <dbl>    <dbl>
    #> 1      1     2  7.14 0.0281      NA        NA       NA
    #> 
    #> $popdrop
    #> # A tibble: 1 x 13
    #>   pat      wt   dof chisq      p f4rank Altai_Neanderth… Vindija.DG feasible
    #>   <chr> <dbl> <int> <dbl>  <dbl>  <int>            <dbl>      <dbl> <lgl>   
    #> 1 00        0     2  7.14 0.0281      1             50.2      -49.2 FALSE   
    #> # … with 4 more variables: best <lgl>, dofdiff <int>, chisqdiff <dbl>,
    #> #   p_nested <dbl>

More documentation
[here](https://uqrmaie1.github.io/admixtools/articles/admixtools.html).

Alternatively, you can launch an interactive browser app on your local
machine:

``` r
run_shiny_admixtools()
```
