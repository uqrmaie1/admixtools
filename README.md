
<!-- README.md is generated from README.Rmd. Please edit that file --->

# admixtools

Lighting fast R implementations of
[AdmixTools](https://github.com/DReichLab/AdmixTools) programs qpAdm,
qpGraph, and more.

## Installation

To install and load `admixtools`, follow these steps:

``` r
devtools::install_github("uqrmaie1/admixtools")
library("admixtools")
```

## Usage

If you already have precomputed f2-statistics, you can fit an admixture
graph like this:

``` r
f2_dir = "/my/f2/directory/"
fit = qpgraph(f2_dir, example_graph)
plot_graph(fit$edges)
```

![example graph](man/figures/graph1.png)

Clearly not a historically accurate model, but it gets the idea across.

<br>

You can also use the f2-statistics to estimate admixture weights:

``` r
target = "Denisova.DG"
left = c("Altai_Neanderthal.DG", "Vindija.DG")
right = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG", "Switzerland_Bichon.SG")
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
    #> 1 Denisova.DG Altai_Neanderthal.DG   49.7  23.3  2.13
    #> 2 Denisova.DG Vindija.DG            -48.7  23.3 -2.08
    #> 
    #> $f4
    #> # A tibble: 36 x 9
    #>    pop1    pop2       pop3   pop4           est      se       z         p weight
    #>    <chr>   <chr>      <chr>  <chr>        <dbl>   <dbl>   <dbl>     <dbl>  <dbl>
    #>  1 Deniso… Altai_Nea… Chimp… Mbuti.DG   0.0101  2.96e-4  34.0   9.79e-253   49.7
    #>  2 Deniso… fit        Chimp… Mbuti.DG   0.00539 7.35e-3   0.733 4.64e-  1   NA  
    #>  3 Deniso… Vindija.DG Chimp… Mbuti.DG   0.0102  2.44e-4  41.7   0.         -48.7
    #>  4 Deniso… Altai_Nea… Chimp… Russia_U…  0.0118  4.53e-4  26.1   5.11e-150   49.7
    #>  5 Deniso… fit        Chimp… Russia_U… -0.00500 1.16e-2  -0.430 6.67e-  1   NA  
    #>  6 Deniso… Vindija.DG Chimp… Russia_U…  0.0122  3.18e-4  38.2   0.         -48.7
    #>  7 Deniso… Altai_Nea… Chimp… Switzerl…  0.0116  4.05e-4  28.7   1.17e-181   49.7
    #>  8 Deniso… fit        Chimp… Switzerl… -0.00430 1.24e-2  -0.346 7.30e-  1   NA  
    #>  9 Deniso… Vindija.DG Chimp… Switzerl…  0.0120  3.00e-4  39.9   0.         -48.7
    #> 10 Deniso… Altai_Nea… Mbuti… Chimp.REF -0.0101  2.96e-4 -34.0   9.79e-253   49.7
    #> # … with 26 more rows
    #> 
    #> $rankdrop
    #> # A tibble: 1 x 7
    #>   f4rank   dof chisq      p dofdiff chisqdiff p_nested
    #> *  <int> <int> <dbl>  <dbl>   <int>     <dbl>    <dbl>
    #> 1      1     2  7.21 0.0272      NA        NA       NA
    #> 
    #> $popdrop
    #> # A tibble: 1 x 13
    #>   pat      wt   dof chisq      p f4rank Altai_Neanderth… Vindija.DG feasible
    #>   <chr> <dbl> <int> <dbl>  <dbl>  <int>            <dbl>      <dbl> <lgl>   
    #> 1 00        0     2  7.21 0.0272      1             49.7      -48.7 FALSE   
    #> # … with 4 more variables: best <lgl>, dofdiff <int>, chisqdiff <dbl>,
    #> #   p_nested <dbl>

More documentation
[here](https://uqrmaie1.github.io/admixtools/articles/admixtools.html).

Alternatively, you can launch an interactive browser app on your local
machine:

``` r
run_shiny_admixtools()
```
