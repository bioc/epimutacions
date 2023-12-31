% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epimutations.R
\name{epimutations}
\alias{epimutations}
\title{Epimutations analysis based on outlier detection methods}
\usage{
epimutations(
  case_samples,
  control_panel,
  method = "manova",
  chr = NULL,
  start = NULL,
  end = NULL,
  epi_params = epi_parameters(),
  maxGap = 1000,
  bump_cutoff = 0.1,
  min_cpg = 3,
  verbose = TRUE
)
}
\arguments{
\item{case_samples}{a GenomicRatioSet object containing the case samples.
See the constructor function \link[minfi]{GenomicRatioSet},
\link[minfi]{makeGenomicRatioSetFromMatrix}.}

\item{control_panel}{a GenomicRatioSet object containing the
control panel (control panel).}

\item{method}{a character string naming the
outlier detection method to be used.
This can be set as: \code{"manova"},
\code{"mlm"}, \code{"iForest"}, \code{"mahdist"},
\code{"quantile"} and \code{"beta"}.
The default is \code{"manova"}.
For more information see \strong{Details}.}

\item{chr}{a character string containing the sequence
names to be analysed. The default value is \code{NULL}.}

\item{start}{an integer specifying the start position.
The default value is \code{NULL}.}

\item{end}{an integer specifying the end position.
The default value is \code{NULL}.}

\item{epi_params}{the parameters for each method.
See the function \link[epimutacions]{epi_parameters}.}

\item{maxGap}{the maximum location gap used in
\link[bumphunter]{bumphunter} method.}

\item{bump_cutoff}{a numeric value of the
estimate of the genomic profile above the
cutoff or below the negative of the
cutoff will be used as candidate regions.}

\item{min_cpg}{an integer specifying the minimum CpGs number in a DMR.}

\item{verbose}{logical. If TRUE additional details about
the procedure will provide to the user.
The default is TRUE.}
}
\value{
The function returns an object of class tibble
containing the outliers regions.
The results are composed by the following columns:
\itemize{
\item \code{epi_id}: systematic name for each epimutation identified.
It provides the name of the used anomaly detection method.
\item \code{sample}: the name of the sample containing the epimutation.
\item \code{chromosome}, \code{start} and \code{end}:
indicate the location of the epimutation.
\item \code{sz}: the window's size of the event.
\item \code{cpg_n}: the number of CpGs in the epimutation.
\item \code{cpg_ids}: the names of CpGs in the epimutation.
\item \code{outlier_score}:
\itemize{
\item For method \code{manova} it provides the approximation
to F-test and the Pillai score, separated by \code{/}.
\item For method \code{mlm} it provides the approximation to
F-test and the R2 of the model, separated by \code{/}.
\item For method \code{iForest} it provides
the magnitude of the outlier score.
\item For method \code{beta} it provides the mean outlier p-value.
\item For methods \code{quantile} and
\code{mahdist} it is filled with NA.
}
\item \code{outlier_direction}: indicates the direction
of the outlier with \code{"hypomethylation"} and \code{"hypermethylation"}
\itemize{
\item For \code{manova}, \code{mlm}, \code{iForest}, and \code{mahdist}
it is computed from the values obtained from bumphunter.
\item For \code{quantile} it is computed from the location
of the sample in the reference distribution (left vs. right outlier).
\item For method \code{beta} it return a NA.
}
\item \code{pvalue}:
\itemize{
\item For methods \code{manova}, \code{mlm}, and \code{iForest}
it provides the p-value obtained from the model.
\item For method \code{quantile}, \code{mahdist} and \code{beta}
is filled with NA.
}
\item \code{adj_pvalue}: for methods with p-value (\code{manova} and
\code{mlm} adjusted p-value with Benjamini-Hochberg based on the total
number of regions detected by Bumphunter.
\item \code{epi_region_id}: Name of the epimutation region as defined
in \code{candRegsGR}.
\item \code{CRE}: cREs (cis-Regulatory Elements) as defined by ENCODE
overlapping the epimutation region. Different cREs are separated by ;.
\item \code{CRE_type}: Type of cREs (cis-Regulatory Elements) as defined
by ENCODE. Different type are separeted by,
and different cREs are separated by ;.
}
}
\description{
The function identifies  differentially methylated regions
in a case sample by comparing it against a control panel.
}
\details{
The function compares a case sample against
a control panel to identify epimutations in the given sample.
First, the DMRs are identified using the
\link[bumphunter]{bumphunter} approach.
After that, CpGs in those DMRs are tested in order to detect regions
with CpGs being outliers.
For that, different outlier detection methods can be selected:
\itemize{
\item Multivariate Analysis of Variance (\code{"manova"}). \link[stats]{manova}
\item Multivariate Linear Model (\code{"mlm"})
\item Isolation Forest (\code{"iForest"}) \link[isotree]{isolation.forest}
\item Robust Mahalanobis Distance (\code{"mahdist"})
\link[robustbase]{covMcd}
\item Quantile distribution (\code{"quantile"})
\item Beta (\code{"beta"})
}

We defined candidate epimutation regions (found in candRegsGR)
based on the 450K array design. As CpGs are not equally distributed
along the genome, only CpGs closer
to other CpGs can form an epimutation.
More information can be found in candRegsGR documentation.
}
\examples{
data(GRset)

#Find epimutations in GSM2562701 sample of GRset dataset

case_samples <- GRset[,11]
control_panel <- GRset[,1:10]
epimutations(case_samples, control_panel, method = "manova")
}
