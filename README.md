amen
----

#### Additive and multiplicative effects network models

[Package website](https://pdhoff.github.io/amen/)

![](circplot.png)

#### About

Additive and multiplicative effects network models (AMEN models) provide
a statistical modeling framework for dyadic network and relational data,
built upon familiar data analysis tools such as linear regression,
random effects models and matrix decompositions. The `amen` package
provides Bayesian model fitting algorithms for AMEN models, and
accommodates a variety of types of network relations, including
continuous, binary and ordinal dyadic variables.

The basic AMEN model is of the form
*y*<sub>*i*, *j*</sub> ∼ *β*<sup>⊤</sup>*x*<sub>*i*, *j*</sub> + *u*<sub>*i*</sub><sup>⊤</sup>*v*<sub>*j*</sub> + *a*<sub>*i*</sub> + *b*<sub>*j*</sub> + *ϵ*<sub>*i*, *j*</sub>
where

-   *y*<sub>*i*, *j*</sub> is the observed dyadic variable being modeled
    and *x*<sub>*i*, *j*</sub> is an observed vector of regressors;

-   *a*<sub>*i*</sub> + *b*<sub>*j*</sub> + *ϵ*<sub>*i*, *j*</sub> is an
    additive random effects term that describes sender and receiver
    variance (such as outdegree and indegree heterogeneity) and dyadic
    correlation;

-   *u*<sub>*i*</sub><sup>⊤</sup>*v*<sub>*j*</sub> is a multiplicative
    random effects term that describes third-order dependence patterns
    (such as transitivity and clustering) and can be estimated and
    analyzed to uncover low-dimensional structure in the network.

#### Installation

``` r
# Current version on GitHub
devtools::install_github("pdhoff/amen") 

# CRAN-approved version on CRAN
install.packages("amen")
```

#### Documentation and citation

A tutorial article and many data analysis examples are available via the
[tutorial](https://github.com/pdhoff/amen/blob/master/inst/tutorial/amen.pdf).
Please cite this as

Hoff, P.D. (2015) “Dyadic data analysis with *amen*”. arXiv:1506.08237.

A review article that provides some mathematical details and derivations
is available on [arXiv](https://arxiv.org/abs/1807.08038). Please cite
this

Hoff, P.D. (2018) “Additive and multiplicative effects network models”.
arXiv:1807.08038.

The first version of the AMEN model appeared in

Hoff, P.D. (2005) “Bilinear mixed-effects models for dyadic data”. JASA
100(469) 286-295.

That version restricted the multiplicative sender and receiver effects
to be equal (*u*<sub>*i*</sub> = *v*<sub>*i*</sub>). The AMEN model in
its current form does not have this restriction. The current AMEN model
first appeared in

Hoff, P.D., Fosdick, B.K., Volfovsky, A. and Stovel, K. (2013)
“Likelihoods for fixed rank nomination networks”. Network Science,
1(3):253–277.

#### Some examples

-   [Modeling a binary outcome](https://pdhoff.github.io/amen/articles/binary_demo.html)

-   [DIY Modeling of a binary outcome](https://pdhoff.github.io/amen/articles/diy_binary_demo.html)

-   [DIY Poisson regression](https://pdhoff.github.io/amen/articles/diy_Poisson_demo.html)
