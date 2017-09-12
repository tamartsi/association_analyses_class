Disease (binary) trait analysis
-------------------------------

-   Disease trait - binary outcome, modeled as *D* = 1 (diseased) or *D* = 0 (disease free).
-   The basic logistic regression model for a binary outcome:
    logit\[*p*(*D*<sub>*i*</sub> = 1)\] = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α*, *i* = 1, …, *n*.
     where here:
    \begin{itemize}
    \item $D_i$ is the disease status of person $i$
    \item  $\mathbf{x}_i$ is a vector of adjusting covariates (age, sex, etc.), $\beta$ is a vector of their effects. 
    \item $g_{i}$ is the dosage or count of the genotype allele of interest.
    \item  $\mbox{logit}(u) = log[u/(1-u)]$, is a function that ensures that estimated disease probabilities - $u$ - will always be in the range $(0,1)$ (while $\text{logit}(u)$ could be anything). 
    \end{itemize}
    Note: there is no "residual". In linear regression, the residual induces the variability. Here, we directly model a probability, which induces a variability.

Disease (binary) trait analysis
-------------------------------

*l**o**g**i**t*\[*p*(*D*<sub>*i*</sub> = 1)\] = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α*, *i* = 1, …, *n*.
\begin{itemize}
\item The basic assumption in this model is that observations are "independent and identically distributed" (i.i.d.). 
\item This does not hold for the HCHS/SOL.
  \begin{itemize}
  \item  So we cannot use the "usual" logistic regression.  
  \item We use mixed models (or GEEs), instead. 
  \end{itemize}
\end{itemize}
Questions:
\begin{enumerate}
\item What will happen if we used logistic regression instead of a logistic mixed model?
\item How can we use logistic regression correctly, assuming we really wanted to?
\end{enumerate}
Disease (binary) trait analysis
-------------------------------

\begin{itemize}
\item The logistic mixed model states that the disease probabilities of people who are somehow close or similar to each other, are more similar to each other than the disease probabilities of people who are not close or similar. 
\item One way to model this is using random effects.
\item For example, if there was one source of such similarities between disease probabilities:
\end{itemize}
*l**o**g**i**t*\[*p*(*D*<sub>*i*</sub> = 1)\] = **x**<sub>*i*</sub><sup>*T*</sup>**β** + *g*<sub>*i*</sub>*α* + *b*<sub>*i*</sub>, *i* = 1, …, *n*,
\begin{itemize}
\item with $b_i$ a random effect, increasing/decreasing the baseline odds of the disease. 
\end{itemize}
Disease (binary) trait analysis
-------------------------------

\begin{itemize}
\item Random effects reflect here  similarity in disease odds across individuals using a correlation structure.
\item As in linear regression, we use matrices to model the correlations between random effects across individuals. 
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
\item Here, the correlation between the random effects of $p_1$ and $p_2$ is 0, and that of $p_1$ and $p_3$ is 0.5. Etc.  
\end{itemize}
Logistic mixed models
---------------------

\begin{itemize}
\item Logistic mixed models are similar to logistic regression, with the addition of random effects.
\item The interpretation is not as "clean" and simple as in linear regression. 
  \begin{itemize}
  \item In linear regression, we used random effects to explicitly model correlation between phenotypes across individuals. 
  \item Here, we do NOT explicitly model the equivalent - correlation between disease probabilities. 
  \end{itemize}
\end{itemize}
Logistic mixed models
---------------------

So how are logistic mixed models practically different from linear mixed models?
\begin{itemize}
\item Variance components are still estimated (but no variance term corresponding to independent errors).
\item There is no straight-forward interpretation of heritability based on variance components. 
\item Computationally, logistic models, and logistic mixed models, take longer to fit (compared to their linear counterparts). 
\item Logistic mixed models for more than a single correlation matrix are implemented in the software GMMAT and R package GENESIS (the same algorithm). 
\end{itemize}
Logistic mixed models
---------------------

\begin{itemize}
\item The GMMAT algorithm uses an approximation, which essentially fits linear mixed models about 4 times, each time for a different "working trait", until both "working traits" and estimated model parameters converge (become about the same as in the previous iteration).  
\item We still fit a null model, as in linear regression.
  \begin{itemize}
  \item It takes about four times longer.
  \end{itemize}
\item We use the null model and the "working traits" to test genotype-disease associations. 
\end{itemize}
Take-home message: the "null model" for binary traits takes 4 times longer to fit than that for quantitative traits. Afterwards computation time is the same.

Logistic vs linear mixed models
-------------------------------

\begin{itemize}
\item In the past, people used linear mixed models instead of logistic mixed models. 
\begin{itemize}
\item Because it saved a lot of computation time. 
\end{itemize}
\item Is it okay to use linear mixed models?
\begin{itemize}
\item Sometimes. But better not! 
\item Basic assumption made by linear mixed models: \textcolor{orange}{residual variance is the same for all people.} 
\item Basic assumption made by logistic model: \textcolor{orange}{if someone has a probability $p$ of disease, the variance of her outcome is $p(1-p)$}.
\end{itemize}
\end{itemize}
Logistic vs linear mixed models
-------------------------------

\begin{itemize}
\item Chen et al. (2016, AJHG) showed in 
\textcolor{purple}{``Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models"} that when 
\begin{itemize}
\item MAF differ between sub-populations in the study
\item Disease prevalence differ between these sub-populations
\end{itemize}
\item LMM test statistics can be either too significant (inflated), or too conservative (deflated). 
\end{itemize}
Logistic vs linear mixed models
-------------------------------

\begin{figure}
\includegraphics[scale = 0.3]{chen_qqplots.pdf} 
\end{figure}
Logistic mixed models - the null model
--------------------------------------

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
outcome <- "disease"
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
              family = "binomial", outcome = outcome, 
              covars = covariates, 
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

    ##      V_HH V_kinship 
    ## 0.1401686 0.0000000

Logistic mixed models - the null model
--------------------------------------

\begin{itemize}
\item Let's look at the results:
\end{itemize}
``` r
nullmod$fixef
```

    ##                     Est         SE        Stat         pval
    ## (Intercept) -13.7912429 1.28407528 115.3521900 6.589170e-27
    ## EV1           0.6019640 0.14250266  17.8441160 2.397595e-05
    ## EV2          -0.4123798 0.15120786   7.4378266 6.386697e-03
    ## sexM          0.2579077 0.27949025   0.8515212 3.561224e-01
    ## age           0.2287610 0.02226677 105.5478697 9.263288e-25
    ## groupuw       0.1138958 0.29014617   0.1540927 6.946546e-01

The logstic mixed model and association testing
-----------------------------------------------

\begin{itemize}
\item After estimating variance components in the ``null model", they are assumed fixed. 
\item We now use this null model object in association testing. 
\end{itemize}
``` r
gds <- GdsGenotypeReader(file.path(dir, 
                           "SISG_snp_dosages.gds"))
# assoc <- assocTestMM(genoData = gds, 
#         nullMMobj = nullmod, family = "binomial")
# try to run! it'll give an error.
```

The logistic mixed model and association testing
------------------------------------------------

``` r
snpAnnot <- getobj(file.path(dir, 
                    "SISG_snp_dosages_snpAnnot.RData"))
genoData <- GenotypeData(gds, 
              snpAnnot=snpAnnot, scanAnnot = scanAnnot)
#assoc <- assocTestMM(genoData = genoData, 
#                       nullMMobj = nullmod)
#                   
# try to run! it'll give an error.
```

The logistic mixed model and association testing
------------------------------------------------

\begin{itemize}
\item We cannot use a Wald test for logistic mixed models 
  \begin{itemize}
  \item Wald test requires estimating genotype effects. In logistic regression, this requires re-estimation of variance components, impossible to do efficiently. 
  \end{itemize}
\item Score tests are "under the null", so they are realistic for GWAS based on logistic mixed models. 
\end{itemize}
``` r
assoc <- assocTestMM(genoData = genoData,
                     nullMMobj = nullmod,
                     test = "Score")
```

    ## Running analysis with 500 Samples and 7463 SNPs

    ## Beginning Calculations...

    ## Block 1 of 2 Completed - 1.513 secs

    ## Block 2 of 2 Completed - 0.5889 secs

The logistic mixed model and association testing
------------------------------------------------

``` r
head(assoc)
```

    ##   snpID chr   n   MAF minor.allele         Score          Var   Score.Stat
    ## 1     1   1 500 0.000            A            NA           NA           NA
    ## 2     2   1 500 0.001            A -9.466802e-05 9.464892e-05 9.468713e-05
    ## 3     3   1 500 0.008            A -6.244369e-01 6.706959e-01 5.813686e-01
    ## 4     4   1 500 0.000            A            NA           NA           NA
    ## 5     5   1 500 0.209            B -4.808782e+00 1.194689e+01 1.935599e+00
    ## 6     6   1 500 0.174            B -6.504122e+00 1.082415e+01 3.908262e+00
    ##   Score.pval
    ## 1         NA
    ## 2 0.99223612
    ## 3 0.44577639
    ## 4         NA
    ## 5 0.16414719
    ## 6 0.04804926

``` r
close(gds)
```

Exercises
---------

\begin{enumerate}
\item Use the results from the GWAS we ran on slide 19. Use the function qqPlot() from the GWASTools package to make a q-q plot figure of the $p$-values from the Score test.
\item Use the following approximation between the Score and Wald test to obtain log(ORs) and ORs for the SNP effects:
  \begin{itemize}
  \item $\beta = \frac{\mbox{score}}{\mbox{var}(\mbox{score})}$
  \item $\mbox{SE}(\beta) = 1/\sqrt{\mbox{var}(\mbox{score})}$
  \item $\mbox{OR} = exp(\beta)$. 
  \end{itemize}
  \begin{itemize}
  \item Which SNP has the highest odds ratio (OR)?
  \end{itemize}
\item Use the function manhattanPlot() from the GWASTools package to generate a Manhattan plots for the Score test $p$-values. Did any of the SNPs achieve genome-wide significance? ($p$-value$=5\times 10{-8}$) array-wide significance?
\end{enumerate}
Exercises
---------

\begin{enumerate}\setcounter{enumi}{3}
\item Which variant is most associated with "disease" among all variants?
\item Run linear mixed model GWAS instead of logistic, treating "disease" as a quantitative trait.
\begin{itemize}
\item ...by first fitting a new null model using fitNullMM()
\item Compare the $p$-values obtained by the two methods. You can use a scatter plot or a q-q plot.
\item Do the two GWAS have the same top SNPs?
\end{itemize}
\item Use the parameter scan.include of function fitNullMM to fit perform association testing only in people from the UNC group.
\begin{itemize}
\item Does the this GWAS have the same top SNP as the GWAS that was run on all participants?
\end{itemize}
\end{enumerate}
