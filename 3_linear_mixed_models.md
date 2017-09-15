Quantitative trait analysis
---------------------------

-   Quantitative trait = continuous outcome. The simplest to analyze.
-   The basic linear regression model for a quantitative outcome:
    *y*<sub>*i*</sub> = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α* + *ϵ*<sub>*i*</sub>, *i* = 1, …, *n*.
     where here:
-   *y*<sub>*i*</sub> is the trait value of person *i*.
-   **x**<sub>*i*</sub> is a vector of adjusting covariates (age, sex, etc.), *β* is a vector of their effects.
-   *g*<sub>*i*</sub> is the genotype dosage or count of the SNP (variant) of interes, *α* its effect.
-   *ϵ*<sub>*i*</sub> is a residual.

Quantitative trait analysis
---------------------------

*y*<sub>*i*</sub> = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α* + *ϵ*<sub>*i*</sub>, *i* = 1, …, *n*.

-   The basic assumption in this linear model is that observations are "independent and identically distributed" (i.i.d.).
-   This does not hold for the HCHS/SOL.
    -   So we cannot use the "usual" linear regression.
    -   We use mixed models (or GEEs), instead.

Questions: 1. What will happen if we used linear regression, i.e. assume, contrary to fact, that participants are i.i.d? 2. How can we use linear regression correctly, if we really wanted to?

Quantitative trait analysis
---------------------------

-   The linear mixed model states that the traits of people who are somehow close or similar to each other, are more similar to each other than the trait values of people who are not close or similar.
    -   I.e. some people's traits are to each other.
-   One way to model these correlations is using random effects.
-   For example, if there was one source of such correlations:

*y*<sub>*i*</sub> = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α* + *b*<sub>*i*</sub> + *ϵ*<sub>*i*</sub>, *i* = 1, …, *n*,
 - with *b*<sub>*i*</sub> a random error - or a random effect - in addition to the i.i.d. errors *ϵ*<sub>*i*</sub>.

Quantitative trait analysis
---------------------------

-   Random effects model the correlation between individuals' trait values.
-   Specifically, one can define a matrix to do that. E.g. a kinship matrix. Or a household matrix!

\begin{eqnarray*}
\text{cor}\left[(b_1, b_2, b_3, \ldots)\right] &=& \begin{array}{cc}
 & \begin{array}{cccc} p_1 & p_2& p_3 & \ldots \end{array} \\
\begin{array}{c}
p_1 \\
p_2\\
p_3\\
\vdots \\
\end{array} & 
  \left(
\begin{array}{cccc}
 1 & 0 & 0.5 & \ldots \\
0 & 1 & 0.5 & \ldots \\
0.5 & 0.5 & 1 & \ldots \\
\vdots & & &  \\
\end{array}
\right) \\
\end{array}
\end{eqnarray*}
-   Here, the correlation between the random effects of persons *p*<sub>1</sub> and *p*<sub>2</sub> is 0, and that of *p*<sub>1</sub> and *p*<sub>3</sub> is 0.5. Etc.

Linear mixed models
-------------------

-   Linear mixed models are similar to linear regression, with the addition of correlation information between peoples' traits.
-   In the HCHS/SOL, we have three correlation matrices: kinship (also called genetic relatedness matrix, GRM), household, and block unit.
-   The kinship matrix was estimated based on the genotyping data (using common variants, MAF≥0.05).
-   The household and block unit matrices were calculated based on who people lived with in the same house, or block unit.

Linear mixed models
-------------------

So how are these correlation matrices used?

-   While the correlation structures (the matrices) are pre-defined, the are not.
-   In linear regression (i.i.d. observations) there is a single residual variance: *σ*<sub>*e*</sub><sup>2</sup>.
    -   It is the variance of the i.i.d residuals: var(*ϵ*<sub>*i*</sub>)=*σ*<sub>*e*</sub><sup>2</sup>.
    -   In other words: a single variance component.
-   In mixed models, there are at least 2 variance components. One for the i.i.d. errors, others correspond to random effects.
    -   In the HCHS/SOL: *σ*<sup>2</sup> = *σ*<sub>*e*</sub><sup>2</sup> + *σ*<sub>*g*</sub><sup>2</sup> + *σ*<sub>*h*</sub><sup>2</sup> + *σ*<sub>*c*</sub><sup>2</sup>

Linear mixed models - variance components
-----------------------------------------

\\begin{itemize} - Variance components are used in two important applications. + Association testing; + Heritability estimation. - Both require estimates of the variance components. - They are estimated by fitting a . + A model that includes the trait, and all adjusting covariates, and the random effects matrices; but not individual genotypes.

Linear mixed models - the null model
------------------------------------

Let's try it!

-   We first load our scanAnnotation object.

``` r
library(GWASTools)
library(GENESIS)
dir <- paste0("/home/postdoc/tsofer/SISG/", 
    "Preparing_simulated_data_2")

scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
scanAnnot
```

    ## An object of class 'ScanAnnotationDataFrame'
    ##   scans: 1 2 ... 500 (500 total)
    ##   varLabels: scanID EV1 ... group (8 total)
    ##   varMetadata: labelDescription

Linear mixed models - the null model
------------------------------------

-   Select outcome, covariates, and load correlation matrices.

``` r
varLabels(scanAnnot)[1:4]
```

    ## [1] "scanID" "EV1"    "EV2"    "sex"

``` r
varLabels(scanAnnot)[5:8]
```

    ## [1] "age"     "trait"   "disease" "group"

``` r
covariates <- c("EV1", "EV2", "sex", "age", "group")
outcome <- "trait"
HH.mat <- getobj(file.path(dir, 
                  "SISG_houshold_matrix.RData"))
kin.mat <- getobj(file.path(dir, 
                  "SISG_relatedness_matrix.RData"))
covMatList <- list(HH = HH.mat, kinship = kin.mat)
```

Linear mixed models - the null model
------------------------------------

``` r
nullmod <- fitNullMM(scanData = scanAnnot,
              outcome = outcome, covars = covariates, 
              covMatList = covMatList, verbose = FALSE)
```

Linear mixed models - the null model
------------------------------------

-   Let's look at the results:

``` r
names(nullmod)
```

    ##  [1] "varComp"           "varCompCov"        "fixef"            
    ##  [4] "betaCov"           "fitted.values"     "resid.marginal"   
    ##  [7] "eta"               "resid.conditional" "logLikR"          
    ## [10] "logLik"            "AIC"               "RSS"              
    ## [13] "workingY"          "model.matrix"      "cholSigmaInv"     
    ## [16] "scanID"            "family"            "converged"        
    ## [19] "zeroFLAG"          "hetResid"

``` r
nullmod$varComp
```

    ##      V_HH V_kinship       V_E 
    ##    0.0000    0.0000  231.7541

Linear mixed models - the null model
------------------------------------

-   Let's look at the results:

``` r
nullmod$fixef
```

    ##                   Est         SE        Stat         pval
    ## (Intercept)  4.213919 2.47521363    2.898325 8.867166e-02
    ## EV1          5.532397 0.68118649   65.962115 4.596743e-16
    ## EV2         -3.191225 0.69292155   21.210297 4.115474e-06
    ## sexM         6.636157 1.37220159   23.388236 1.323857e-06
    ## age          3.771601 0.04967911 5763.733606 0.000000e+00
    ## groupuw     -4.847355 1.39223675   12.122256 4.982359e-04

Linear mixed models - the null model
------------------------------------

-   Our simulated trait "trait" unfortunately doesn't seem to be very heritable.
-   Let's simulate another outcomes to make it more interesting...

``` r
require(mvtnorm)
n <- nrow(kin.mat)
new.trait <- 
  nullmod$model.matrix %*% matrix(c(4, 5, -2, 1, 4,2)) + 
    matrix(rmvnorm(n = 1, mean = rep(0, n), 
            sigma = diag(rep(120, n)) + 
                      80*kin.mat + 40*HH.mat))
scanAnnot$new.trait <- as.numeric(new.trait)
```

Linear mixed models - the null model
------------------------------------

``` r
set.seed(101)
nullmod <- fitNullMM(scanData = scanAnnot,
              outcome = "new.trait", 
              covars = covariates, 
              covMatList = covMatList, verbose = FALSE)

nullmod$varComp
```

    ##      V_HH V_kinship       V_E 
    ## 125.81207  86.84574  47.53430

``` r
varCompCI(nullmod, prop = TRUE)
```

    ##           Proportion    Lower 95  Upper 95
    ## V_HH       0.4835353  0.04957174 0.9174989
    ## V_kinship  0.3337755 -0.72562827 1.3931792
    ## V_E        0.1826892 -0.93299124 1.2983697

The linear mixed model and heritability
---------------------------------------

-   The proportion of variance due to kinship/genetic relatedness is .
    -   AKA narrow-sense heritability.
    -   The heritability of "trait" is estimated to be 33%, with 95% confidence interval (-73, 139)%.
    -   To test heritability we can use the confidence intervals - if they are calculated correctly(!), or the likelihood ratio test.
-   The simulated data set has 500 people, which is very small.
-   Therefore, variance components are not well estimated,
-   and the confidence interval of the heritabability includes impossible values.
    -   Negative values, and larger than 100...
    -   There are methods to calculate feasible CIs.

The linear mixed model and association testing
----------------------------------------------

-   After estimating variance components in the \`\`null model", they are assumed fixed.
-   We now use this null model object in association testing.
    -   Note: it can take a long time to estimate variance components, so doing it once (instead of separately for every genetic variant) saves a lot of time.

``` r
gds <- GdsGenotypeReader(file.path(dir, 
                           "SISG_snp_dosages.gds"))
# assoc <- assocTestMM(genoData = gds, 
#                   nullMMobj = nullmod)
# try to run! it'll give an error.
```

The linear mixed model and association testing
----------------------------------------------

``` r
snpAnnot <- getobj(file.path(dir, 
                    "SISG_snp_dosages_snpAnnot.RData"))
genoData <- GenotypeData(gds, 
              snpAnnot=snpAnnot, scanAnnot = scanAnnot)
assoc <- assocTestMM(genoData = genoData, 
                     nullMMobj = nullmod)
```

    ## Running analysis with 500 Samples and 7463 SNPs

    ## Beginning Calculations...

    ## Block 1 of 2 Completed - 1.261 secs

    ## Block 2 of 2 Completed - 0.9561 secs

The linear mixed model and association testing
----------------------------------------------

``` r
head(assoc)
```

    ##   snpID chr   n   MAF minor.allele        Est        SE  Wald.Stat
    ## 1     1   1 500 0.000            A         NA        NA         NA
    ## 2     2   1 500 0.001            A -2.3574390 16.248231 0.02105081
    ## 3     3   1 500 0.008            A -0.6907384  5.655461 0.01491733
    ## 4     4   1 500 0.000            A         NA        NA         NA
    ## 5     5   1 500 0.209            B  1.4273077  1.211778 1.38736071
    ## 6     6   1 500 0.174            B  0.7567469  1.304453 0.33654615
    ##   Wald.pval
    ## 1        NA
    ## 2 0.8846406
    ## 3 0.9027909
    ## 4        NA
    ## 5 0.2388513
    ## 6 0.5618297

``` r
close(gds)
```

Inflation in Genome-Wide Association Studies
--------------------------------------------

-   Fundamental assumption in GWAS:
    -   Most genetic variants are not associated with the outcome.
-   Test statistics are mostly distributed \`\`under the null"
-   $\_{gc} = $
    -   Ideally, *λ* = 1.
    -   Because the bulk of the association are null.
-   q-q plots are used to evaluate inflation. ![qq plot example](/figures/qq_fig.png)

Exercises
---------

1.  Use the results from the GWAS that we ran on slide 18, and the function qqPlot() from the GWASTools package to make a q-q plot of the *p*-values from the Wald test.
2.  What is the inflation factor? use the R code to obtain the median of the expected distribution of the test statistics, and $ to obtain the median of the observed test statistics.
3.  Can you evaluate whether the GWAS is too inflated or deflated?
4.  Use the function manhattanPlot() from the GWASTools package to make a Manhattan plots for these *p*-values.

Exercises
---------

1.  Set the significance threshold line in the Manhattan plot to be 0.05 divided by the number of tested variants (Bonferroni correction).
2.  Show results in the Manhattan plot only for variants with imputation quality ("info") at least 0.8, or genotyped.
3.  Which variant is most associated with "trait" among all variants (according to *p*-value)?
4.  Use the parameter snp.include of function assocTestMM() to test only variants in positions 1029889 - 2136826 on chromosome 1.
    -   Which variant has the most significant *p*-value?

5.  Use the parameter scan.include of function fitNullMM() to perform association testing only in people from the UW group.
    -   Is the most significant variant the same as before?
