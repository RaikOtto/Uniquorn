\name{NEWS}
\title{NEWS file for the \pkg{Uniquorn} bioconductor package }
\section{Changes in version 2.1.3 [2018-10-03]}{
    \subsection{Bioconductor compliance}{
        \itemize{
            \item Minor bugfixes regarding the removal of CCLs and reference libraries
        }
    }
}
\section{Changes in version 1.99.3 [2018-01-03]}{
    \subsection{Bioconductor compliance}{
        \itemize{
            \item Minor modifications to comply with Bioconductor regulations
        }
    }
}
\section{Changes in version 1.99.2 [2018-20-02]}{
    \subsection{Added RNA-seq and panel-seq capability}{
        \itemize{
            \item Some more minor bugfixes
        }
    }
}
\section{Changes in version 1.99.1 [2018-12-02]}{
    \subsection{Added RNA-seq and panel-seq capability}{
        \itemize{
            \item Minor bugfixes
        }
    }
}
\section{Changes in version 1.99.0 [2018-15-01]}{
    \subsection{Added RNA-seq and panel-seq capability}{
        \itemize{
            \item Modified the package to work with large scale RNA-seq and small scale panel-seq data
            \item Please note that this might break custom libraries. In case of questions please contact the author
        }
    }
}
\section{Changes in version 1.4.2 [2017-15-05]}{
    \subsection{Added critical scripts}{
        \itemize{
            \item Added a script to parse vcf files in a non-parallelized fashion.
        }
    }
}
\section{Changes in version 1.4.1 [2017-12-05]}{
    \subsection{Fixed yet another penalty function bug}{
        \itemize{
            \item Fixed critical bug which disabled identifying CCLs in single-thread mode.
        }
    }
}
\section{Changes in version 1.3.1 [2016-28-11]}{
    \subsection{Fixed yet another penalty function bug}{
        \itemize{
            \item Fixed rather major bug that led to too high penalty scores.
        }
    }
}
\section{Changes in version 1.1.3 [2016-28-11]}{
    \subsection{Fixed penalty function bug}{
        \itemize{
            \item Fixed minor bug that led to incorrect penalty score calculation.
        }
    }
}
\section{Changes in version 1.1.2 [2016-28-11]}{
    \subsection{Change of penalty function}{
        \itemize{
            \item The function that measures random background-noise matches and calculates a threshold above the background strength has been updated to be more strict.
        }
    }
}
\section{Changes in version 1.1.1 [2016-10-18]}{
    \subsection{Change of default inclusion threshold}{
        \itemize{
            \item the new default inclusion weight has been set to 0.5 to optimally fit the new probabilistic model of version 1.0.9.
        }
    }
}
\section{Changes in version 1.0.9 [2016-10-11]}{
    \subsection{Re-engineering of probabilistic model}{
        \itemize{
            \item The likelihood to observe a mutation/ variant in a query that matches a reference CCL's variant has been changed. The likelihood was 1 / #all variants in library and is now set to # number of variants in reference CCL / #all variants in library
        }
    }
}
\section{Changes in version 1.0.8 [2016-06-17]}{
    \subsection{Optimization confidence score}{
        \itemize{
            \item Fixed the installation of the package
        }
    }
}
\section{Changes in version 1.0.7 [2016-06-13]}{
    \subsection{Optimization confidence score}{
        \itemize{
            \item Optimized the way confidence scores are calculated to increase the sensitivity and specificity of the method 
        }
    }
}
\section{Changes in version 1.0.6 [2016-06-05]}{
    \subsection{Optimization of default confidence score}{
        \itemize{
            \item Adjusted default confidence score to optimal threshold 
        }
    }
}
\section{Changes in version 1.0.5 [2016-05-30]}{
    \subsection{Minor Bugfixes}{
        \itemize{
            \item Fixed errors connected to adding CLs
        }
    }
}
\section{Changes in version 1.0.4 [2016-05-27]}{
    \subsection{Introduction of confidence score}{
        \itemize{
            \item Confidence score is the central threshold now
            \item Default is 3.0
            \item The score is the negative log e of the q-value
        }
    }
}
\section{Changes in version 1.0.3 [2016-05-18]}{
    \subsection{Changes to treshold calculation}{
        \itemize{
            \item Fixed error in threshold calculation
        }
    }
}
\section{Changes in version 1.0.1 [2016-05-09]}{
    \subsection{Minor update of BED files}{
        \itemize{
            \item Extended length of BED file mutations shown in the IGV
        }
    }
}
\section{Changes in version 0.99.16 [2016-05-04]}{
    \subsection{Further minor update of statistics module}{
        \itemize{
            \item Further minor fixes to statistics
        }
    }
}
\section{Changes in version 0.99.15 [2016-05-03]}{
    \subsection{Minor Update of statistics}{
        \itemize{
            \item Minor fixes to statistics
        }
    }
}
\section{Changes in version 0.99.14 [2016-05-03]}{
    \subsection{Major Update of statistics}{
        \itemize{
            \item Replaced the absolute and relative cutoff with a p an q-value based on a binomial distribution
        }
    }
}
\section{Changes in version 0.99.13 [2016-04-08]}{
    \subsection{Fixed unit tests}{
        \itemize{
            \item Further update test suite
        }
    }
}
\section{Changes in version 0.99.12 [2016-04-08]}{
    \subsection{Fixed unit tests}{
        \itemize{
            \item Update testing suite
        }
    }
}
\section{Changes in version 0.99.11 [2016-04-01]}{
    \subsection{Fixed unit tests}{
        \itemize{
            \item Fixed error in unit tests
        }
    }
}
\section{Changes in version 0.99.10 [2016-03-30]}{
    \subsection{Added Unit tests}{
        \itemize{
            \item Added Unit test feature
        }
    }
}
\section{Changes in version 0.99.9 [2016-03-29]}{
    \subsection{Improved documentation}{
        \itemize{
            \item Minor documentation bugfixes.
        }
    }
}
\section{Changes in version 0.99.8 [2016-03-25]}{
    \subsection{Fixed Bug}{
        \itemize{
            \item General bugfixing.
        }
    }
}
\section{Changes in version 0.99.7 [2016-03-24]}{
    \subsection{Fixed Bug}{
        \itemize{
            \item Minor bugfix concerning output of CL identifications.
        }
    }
}
\section{Changes in version 0.99.6 [2016-03-23]}{
    \subsection{Fixed Bug}{
        \itemize{
            \item Minor bugfix concerning data parsing of CellMiner data.
        }
    }
}
\section{Changes in version 0.99.5 [2016-03-22]}{
    \subsection{Fixed Bug}{
        \itemize{
            \item Fixed a bug that lead to incorrect calculation of CoSMIC CLP CL mutation weights.
        }
    }
}
\section{Changes in version 0.99.4 [2016-03-16]}{
    \subsection{Added features}{
        \itemize{
            \item Minimum amount of matching mutations required for a positive prediction between query and training sample can now be manually adjusted in the 'identify_vcf_file' function using the parameter 'miminum_matching_mutations'
        }
    }
}

\section{Changes in version 0.99.3 [2016-03-02]}{
    \subsection{Fixed bugs}{
        \itemize{
            \item Fixed vigniette
        }
    }
}

\section{Changes in version 0.99.2 [2016-03-02]}{
    \subsection{Fixed bugs}{
        \itemize{
            \item Moved NEWS.Rd from 'Uniquorn' to 'Uniquorn/inst'
        }
    }
}

\section{Changes in version 0.99.1 [2016-03-02]}{
    \subsection{Fixed bugs}{
        \itemize{
            \item Fixed NEWS section
            \item Fixed unit testing
            \item Fixed the Vignette: replaced the github-dependent installation with the biocLite installation
        }
    }
}

\section{Changes in version 0.99.0 [2015-01-27]}{
    \subsection{New features}{
        \itemize{
            \item Identify cancer cell lines
            \item Show which cancer cell line are contained
            \item Show which mutations are annotated for selected cancer cell lines
            \item Show which mutaitons are overall included
            \item Parse cancer cell line custom data -> add your own samples and identify these
        }
    }
}
