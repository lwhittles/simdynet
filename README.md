
<!-- README.md is generated from README.Rmd. Please edit that file -->
simplasex
=========

The goal of simplasex is to simulate realistic dynamic sexual networks.

Installation
------------

You can install simplasex from github with:

``` r
devtools::install_github("lwhittles/simplasex")
```

The package can then be loaded using:

``` r
library(simplasex)
```

Example
-------

This is a basic example of usage.

``` r
test_og <- sim_dynamic_sn(N = 1e3, gamma = 1.8, k0 = 0.5, phi = 1e4,
                          t = 1+1/365, max.iter = 1e6)
#> Time difference of 0.573673 secs
```

More information and getting help
---------------------------------

For more detailed examples of how to use simplasex, see the vignettes [here](https://github.com/xavierdidelot/simplasex/tree/master/vignettes). See also the help included in the package using the R command `help(package='simplasex')`.

If you have any problem or question please create an issue [here](https://github.com/lwhittles/simplasex/issues) or get in touch by emailing `l.whittles14@imperial.ac.uk` or `xavier.didelot@gmail.com`
