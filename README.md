
<!-- README.md is generated from README.Rmd. Please edit that file --->

# ADMIXTOOLS 2.0

A new, lightning fast implementation of
[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools).

## Overview

ADMIXTOOLS is a set of programs which use f-statistics to infer
demographic models. It has been used in countless publications to test
whether populations form clades (*qpDstat*, *qpWave*), to estimate
ancestry proportions (*qpAdm*), and to fit admixture graphs (*qpGraph*).

ADMIXTOOLS 2.0 provides the same functionality in a new look, and it’s
orders of magnitude faster. This is achieved mostly through separating
the computation of f2-statistics from all other computations. In the
example below, rendering the plot takes much longer than computing the
fit of a new *qpGraph* model:

![app demo](man/figures/shinyapp1.gif)

## Features

  - Much faster than the original ADMIXTOOLS software
  - Simple R command line interface
  - Even simpler point-and-click interface
  - Several new features and methodological innovations that make it
    easier to find robust models:
      - Automated and semi-automated admixture graph inference
      - Simultaneous exploration of hundreds of *qpAdm* models
      - Unbiased comparison of any two *qpGraph* models using
        out-of-sample scores
      - Jackknife and bootstrap confidence intervals for any *qpAdm*,
        *qpWave*, and *qpGraph* parameters
      - Detailed output for each fitted model
  - Full support for genotype data in (PACKED)ANCESTRYMAP/EIGENSTRAT
    format and PLINK format
  - Wrapper functions around the original ADMIXTOOLS software (see also
    [admixr](https://bodkan.net/admixr/index.html))
    <!-- * Simple interface with [msprime](https://msprime.readthedocs.io/en/stable/index.html) for simulating under a given admixture graph -->
  - [Extensive
    documentation](https://uqrmaie1.github.io/admixtools/articles/admixtools.html)
  - New features available [on
    request](mailto:rmaier@broadinstitute.org)\!

## Installation

To install and load ADMIXTOOLS 2.0, start R (version 3.5 or higher) and
follow these steps:

``` r
install.packages("devtools") # if "devtools" is not installed already
devtools::install_github("uqrmaie1/admixtools")
library("admixtools")
```

The above commands will install all R package dependencies which are
required to run ADMIXTOOLS 2.0 on the command line. For the interactive
app, additional packages are required, which can be installed like this:

``` r
devtools::install_github("uqrmaie1/admixtools", dependencies = TRUE)
```

If you encounter any problems during the installation, this is most
likely because some of the required R packages cannot be installed. If
that happens, try re-installing some of the larger packages, and pay
attention to any error message:

``` r
install.packages("Rcpp")
install.packages("tidyverse")
install.packages("igraph")
install.packages("plotly")
```

If you get the following error on Linux `Error: package or namespace
load failed for 'admixtools' in dyn.load(file, DLLpath = DLLpath, ...):`
try adding the following two lines to the file `~/.R/Makevars` (and
create the file first if it doesn’t exist)

    #LAPACK_LIBS=-llapack
    PKG_LIBS = $(LAPACK_LIBS)

If the installation still fails, please [contact
me](mailto:rmaier@broadinstitute.org).

## Usage

First we need to extract f2-statistics from genotype files. These
f2-statistics will be written to disk so that the slow part of the
computation does not have to be repeated.

``` r
genotype_data = "/my/geno/prefix"
f2_dir = "/my/f2/directory/"
extract_f2(genotype_data, f2_dir)
```

After that, we can fit an admixture graph like this:

``` r
f2_data = f2_from_precomp(f2_dir)
fit = qpgraph(f2_data, example_graph)
plot_graph(fit$edges)
```

![example graph](man/figures/graph1.png)

Clearly not a historically accurate model, but it gets the idea across.

<br>

We can also use the f2-statistics to estimate admixture weights:

``` r
target = "Denisova.DG"
left = c("Altai_Neanderthal.DG", "Vindija.DG")
right = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG", "Switzerland_Bichon.SG")
```

``` r
qpadm(f2_data, target, left, right)$weights
```

    #> # A tibble: 2 x 5
    #>   target      left                 weight    se     z
    #>   <chr>       <chr>                 <dbl> <dbl> <dbl>
    #> 1 Denisova.DG Altai_Neanderthal.DG   49.6  23.3  2.13
    #> 2 Denisova.DG Vindija.DG            -48.6  23.3 -2.08

<br>

Or we can get f4-statistics from the f2-statistics:

``` r
f4(f2_data)
```

    #> # A tibble: 105 x 8
    #>    pop1        pop2         pop3    pop4            est      se      z         p
    #>    <chr>       <chr>        <chr>   <chr>         <dbl>   <dbl>  <dbl>     <dbl>
    #>  1 Altai_Nean… Chimp.REF    Deniso… Mbuti.DG     0.0196 6.07e-4   32.4 1.32e-229
    #>  2 Altai_Nean… Denisova.DG  Chimp.… Mbuti.DG    -0.0129 3.64e-4  -35.6 2.22e-277
    #>  3 Altai_Nean… Mbuti.DG     Chimp.… Denisova.DG -0.0326 5.22e-4  -62.5 0.       
    #>  4 Altai_Nean… Chimp.REF    Deniso… Russia_Ust…  0.0180 6.87e-4   26.3 6.43e-152
    #>  5 Altai_Nean… Denisova.DG  Chimp.… Russia_Ust… -0.0152 4.46e-4  -34.0 4.67e-254
    #>  6 Altai_Nean… Russia_Ust_… Chimp.… Denisova.DG -0.0332 5.55e-4  -60.0 0.       
    #>  7 Altai_Nean… Chimp.REF    Deniso… Switzerlan…  0.0181 6.63e-4   27.3 1.09e-164
    #>  8 Altai_Nean… Denisova.DG  Chimp.… Switzerlan… -0.0150 4.64e-4  -32.3 6.06e-229
    #>  9 Altai_Nean… Switzerland… Chimp.… Denisova.DG -0.0331 5.74e-4  -57.7 0.       
    #> 10 Altai_Nean… Chimp.REF    Deniso… Vindija.DG  -0.0771 6.98e-4 -110.  0.       
    #> # … with 95 more rows

## Interactive browser app

ADMIXTOOLS 2.0 also has a simple point-and-click interface. This makes
it easy to explore many *qpAdm* or *qpGraph* models at the same time,
for example by allowing you to build and change admixture graphs
interactively. Typing the following command in the R console launches
the app:

``` r
run_shiny_admixtools()
```

![app demo](man/figures/shinyapp2.gif)

## Documentation

One of the design goals behind ADMIXTOOLS 2.0 is to make the algorithms
more transparent, so that the path from genotype data to conclusions
about demographic history is easier to follow.

To this end, all programs and parameters are (or should be) explained
[here](https://uqrmaie1.github.io/admixtools/articles/admixtools.html)

In addition to that, many of the core functions are implemented twice:
Once in C++ for performance (used by default), and another time in R
where it is easier to trace the computations step by step.

## Contact

For questions, feature requests, and bug reports, please contact Robert
Maier <rmaier@broadinstitute.org>.

## See also

  - [ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) The original
    ADMIXTOOLS software
  - [admixr](https://bodkan.net/admixr/index.html) An R package with
    ADMIXTOOLS wrapper functions and many useful tutorials
  - [admixturegraph](https://github.com/mailund/admixture_graph) An R
    package for automatic graph inference
  - [miqoGraph](https://github.com/juliayyan/PhylogeneticTrees.jl) A
    Julia package for automatic graph inference
  - [MixMapper](http://cb.csail.mit.edu/cb/mixmapper/) Another method to
    infer admixture graphs
  - [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)
    Another method to infer admixture
    graphs
  - [Legofit](http://content.csbs.utah.edu/~rogers/src/legofit/index.html)
    A program to estimate the history of population size, subdivision,
    and gene flow
