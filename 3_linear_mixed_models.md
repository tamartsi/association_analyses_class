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
\begin{itemize}
\item The basic assumption in this linear model is that observations are "independent and identically distributed" (i.i.d.). 
\item This does not hold for the HCHS/SOL.
\begin{itemize}
\item  So we cannot use the "usual" linear regression.  
\item We use mixed models (or GEEs), instead. 
\end{itemize}
\end{itemize}
Questions:
\begin{enumerate}
\item What will happen if we used linear regression, i.e. assume, contrary to fact, that participants are i.i.d?
\item How can we use linear regression correctly, if we really wanted to?
\end{enumerate}
Quantitative trait analysis
---------------------------

\begin{itemize}
\item The linear mixed model states that the traits of people who are somehow close or similar to each other, are more similar to each other than the trait values of people who are not close or similar. 
\begin{itemize}
    \item I.e. some people's traits are \textcolor{purple}{correlated} to each other. 
    \end{itemize}
\item One way to model these correlations is using random effects.
\item For example, if there was one source of such correlations:
\end{itemize}
*y*<sub>*i*</sub> = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α* + *b*<sub>*i*</sub> + *ϵ*<sub>*i*</sub>, *i* = 1, …, *n*,
\begin{itemize}
\item with $b_i$ a random error - or a random effect - in addition to the i.i.d. errors $\epsilon_i$.  
\end{itemize}
Quantitative trait analysis
---------------------------

\begin{itemize}
\item Random effects model the correlation between individuals' trait values.
\item Specifically, one can define a matrix to do that. E.g. a kinship matrix. Or a household matrix!
\end{itemize}
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
\end{eqnarray*}\begin{itemize}
\item Here, the correlation between the random effects of persons $p_1$ and $p_2$ is 0, and that of $p_1$ and $p_3$ is 0.5. Etc.  
\end{itemize}
Linear mixed models
-------------------

\begin{itemize}
\item Linear mixed models are similar to linear regression, with the addition of correlation information between peoples' traits. 
\item In the HCHS/SOL, we have three correlation matrices: kinship (also called genetic relatedness matrix, GRM), household, and block unit. 
\item The kinship matrix was estimated based on the genotyping data (using common variants, MAF$\geq 0.05$).
\item The household and block unit matrices were calculated based on who people lived with in the same house, or block unit. 
\end{itemize}
Linear mixed models
-------------------

So how are these correlation matrices used?
\begin{itemize}
\item While the correlation structures (the matrices) are pre-defined, the \textcolor{purple}{ variance components} are not. 
\item In linear regression (i.i.d. observations) there is a single residual variance: $\sigma_e^2$. 
  \begin{itemize}
  \item It is the variance of the i.i.d residuals: $\mbox{var}(\epsilon_i) = \sigma_e^2$. 
  \item In other words: a single variance component. 
  \end{itemize}
\item In mixed models, there are at least 2 variance components. One for the i.i.d. errors, others correspond to random effects. 
  \begin{itemize}
  \item  In the HCHS/SOL: $\sigma^2 = \sigma^2_e + \sigma^2_g + \sigma^2_h + \sigma^2_c$
  \end{itemize}
\end{itemize}
Linear mixed models - variance components
-----------------------------------------

\begin{itemize}
\item Variance components are used in two important applications. 
  \begin{itemize}
  \item Association testing;
  \item Heritability estimation.
  \end{itemize}
\item Both require estimates of the variance components. 
\item They are estimated by fitting a \textcolor{purple}{null model}.
  \begin{itemize}
  \item A model that includes the trait, and all adjusting covariates, and the random effects matrices; but not individual genotypes. 
  \end{itemize}
\end{itemize}
Linear mixed models - the null model
------------------------------------

Let's try it!
\begin{itemize}
\item We first load our scanAnnotation object. 
\end{itemize}
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

\begin{itemize}
\item Select outcome, covariates, and load correlation matrices.
\end{itemize}
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

\begin{itemize}
\item Let's look at the results:
\end{itemize}
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

\begin{itemize}
\item Let's look at the results:
\end{itemize}
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

\begin{itemize}
\item Our simulated trait "trait" unfortunately doesn't seem to be very heritable.
\item Let's simulate another outcomes to make it more interesting... 
\end{itemize}
``` r
require(mvtnorm)
```

    ## Loading required package: mvtnorm

``` r
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
    ##   0.00000  21.72032 210.16076

``` r
varCompCI(nullmod, prop = TRUE)
```

    ##           Proportion   Lower 95 Upper 95
    ## V_HH      0.00000000         NA       NA
    ## V_kinship 0.09367009 -1.0453772 1.232717
    ## V_E       0.90632991 -0.2327174 2.045377

The linear mixed model and heritability
---------------------------------------

\begin{itemize}
\item The proportion of variance due to kinship/genetic relatedness is \textcolor{purple}{heritability}.
  \begin{itemize}
  \item AKA narrow-sense heritability.
  \item The heritability of "trait" is estimated to be 9\%, with 95\% confidence interval (-105, 123)\%.
  \item To test heritability we can use the confidence intervals - if they are calculated correctly(!), or the likelihood ratio test. 
  \end{itemize}
\item The simulated data set has 500 people, which is very small. 
\item Therefore, variance components are not well estimated,
\item and the confidence interval of the heritabability includes impossible values.
  \begin{itemize}
  \item Negative values, and larger than 100...
  \item There are methods to calculate feasible CIs. 
  \end{itemize}
\end{itemize}
The linear mixed model and association testing
----------------------------------------------

\begin{itemize}
\item After estimating variance components in the ``null model", they are assumed fixed. 
\item We now use this null model object in association testing. 
  \begin{itemize}
  \item Note: it can take a long time to estimate variance components, so doing it once (instead of separately for every genetic variant) saves a lot of time. 
  \end{itemize}
\end{itemize}
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

    ## Block 1 of 2 Completed - 2.234 secs

    ## Block 2 of 2 Completed - 0.9063 secs

The linear mixed model and association testing
----------------------------------------------

``` r
head(assoc)
```

    ##   snpID chr   n   MAF minor.allele          Est        SE   Wald.Stat
    ## 1     1   1 500 0.000            A           NA        NA          NA
    ## 2     2   1 500 0.001            A -20.89367957 15.320739 1.859818001
    ## 3     3   1 500 0.008            A  -3.02127924  5.452748 0.307008870
    ## 4     4   1 500 0.000            A           NA        NA          NA
    ## 5     5   1 500 0.209            B  -0.28969002  1.166012 0.061724953
    ## 6     6   1 500 0.174            B   0.07804366  1.244438 0.003933045
    ##   Wald.pval
    ## 1        NA
    ## 2 0.1726458
    ## 3 0.5795215
    ## 4        NA
    ## 5 0.8037901
    ## 6 0.9499943

``` r
close(gds)
```

Inflation in Genome-Wide Association Studies
--------------------------------------------

\begin{itemize}
\item Fundamental assumption in GWAS:
    \begin{itemize}
    \item Most genetic variants are not associated with the outcome.
    \end{itemize}
\item Test statistics are mostly distributed ``under the null"
\item $\lambda_{gc} = \frac{\mbox{median}(\mbox{observed test statistics})}{\mbox{median}(\mbox{expected distribution of test statistics})} $
\begin{itemize}
\item Ideally, $\lambda = 1$.
\item Because the bulk of the association are null. 
\end{itemize}
\item q-q plots are used to evaluate inflation.
\end{itemize}\begin{center}
\includegraphics[scale = 0.3]{qq_fig.png}
\end{center}
Exercises
---------

\begin{enumerate}
\item Use the results from the GWAS that we ran on slide 18, and the function qqPlot() from the GWASTools package to make a q-q plot of the $p$-values from the Wald test. 
\item What is the inflation factor? use the R code \texttt{pchisq(0.5, df = 1, lower.tail = FALSE)} to obtain the median of the expected distribution of the test statistics, and \texttt{median(assoc}\$\texttt{Wald.Stat, na.rm = TRUE)} to obtain the median of the observed test statistics.
\item Can you evaluate whether the GWAS is too inflated or deflated?
\item Use the function manhattanPlot() from the GWASTools package to make a Manhattan plots for these $p$-values. 
\end{enumerate}
Exercises
---------

\begin{enumerate}\setcounter{enumi}{4}
\item Set the significance threshold line in the Manhattan plot to be 0.05 divided by the number of tested variants (Bonferroni correction).
\item Show results in the Manhattan plot only for variants with imputation quality ("info") at least 0.8, or genotyped.
\item Which variant is most associated with "trait" among all variants (according to $p$-value)?
\item Use the parameter snp.include of function assocTestMM() to test only variants in positions 1029889 - 2136826 on chromosome 1. 
    \begin{itemize}
    \item Which variant has the most significant $p$-value?
    \end{itemize}
\item Use the parameter scan.include of function fitNullMM() to perform association testing only in people from the UW group. 
    \begin{itemize}
    \item Is the most significant variant the same as before?
    \end{itemize}
\end{enumerate}
