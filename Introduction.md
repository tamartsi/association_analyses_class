Overview
--------

-   In this module we will present and discuss methods and softwares for genetic association studies.
-   We will focus on studies with complex population structure.
-   We will follow the Hispanic Community Health Study/Study of Latinos (HCHS/SOL) as a case study.

Overview
--------

We will...

1.  Present the HCHS/SOL.
2.  Present a simulated data set that mimics the HCHS/SOL (simpler).
3.  Learn how to run
    -   Linear mixed model GWAS.
    -   Logistic mixed model GWAS.
    -   GEE based GWASs.

4.  Discuss and practice possible ways to deal with heterogeneity within the study.
    -   In both the mixed models and the GEE frameworks.

5.  Learn to preform generalization analysis.
6.  Run admixture mapping GWAS.
7.  Discuss and practice GxE analyses.
    -   In both the mixed models and the GEE frameworks.

See module schedule for details!

Resources
---------

The class site contains

-   PDF slide sets.
-   Data sets.
-   Instructions for software installation.
-   Published manuscripts of interest.

Instruction team: Tamar Sofer
-----------------------------

\begin{figure}
\includegraphics[scale=0.5]{tsofer_pic.png}
\end{figure}
-   Research scientist at the University of Washington, Department of Biostatistics.
    -   Future position at Harvard Medical School/Brigham and Women's Hospital.
-   Worked at the Genetic Analysis Center of the HCHS/SOL.

Instruction team: Kari North
----------------------------

\begin{figure}
\includegraphics[scale=0.2]{knorth_pic.png}
\end{figure}
-   Professor of Epidemiology at UNC.
-   Leads the Population Architecture using Genomics and Epidemiology (PAGE) study.

Instruction team: Mariaelisa (Misa) Graff
-----------------------------------------

\begin{figure}
\includegraphics[scale=0.2]{misa_pic.png}
\end{figure}
-   Research Assistant Professor of Epidemiology at UNC.
-   Worked on many genetic epidemiology studies.

Now let's download and install softwares!
-----------------------------------------

-   First, download and install the latest R version.
-   Install Bioconductor R packages:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("GWASTools")
biocLite("GENESIS")
biocLite("gdsfmt")
```

-   If you cannot install the latest R version, you may have to install GENESIS package manually to have all functionalities!

Now let's download and install softwares!
-----------------------------------------

-   Install R packages from GitHub:

``` r
install.packages("devtools")
library("devtools")
install_github("tamartsi/generalize@Package_update", 
         subdir = "generalize")
install_github("tamartsi/MetaCor")
```

Now let's download and install softwares!
-----------------------------------------

-   Other R packages.

``` r
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("RColorBrewer")
```

Download data sets
------------------

-   Download data sets from the module website .
-   And save all data sets in the same folder.
-   When using R code, we will always start by setting our working directory using a command that looks like:

``` r
dir <- "/mycomputer/variousfolders/module12_folder"
# put your own!
```

-   Save this command somewhere, with YOUR working directory (the path to the folder you use) for an easy use in the next sessions.
