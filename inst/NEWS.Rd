\name{NEWS}
\title{News for \R Package \pkg{rsae}}
\encoding{UTF-8}
\section{CHANGES in rsae VERSION 0.1-6 (2021-06-20)}{
    \subsection{BUG FIXES}{
        \itemize{
            \item fixed namespace issues
            \item fixed DESCRIPTION file
            \item replaced 'require' in function '.initmethod ' with
                'requireNamespace'
        }
    }
    \subsection{MISC}{
        \itemize{
            \item updated CITATION file
            \item cleaned documentation
            \item combined several R source files into one source file
            \item source repository moved from R-Forge to GitHub
        }
    }
}
\section{CHANGES in rsae VERSION 0.1-5 (2014-02-13)}{
    \subsection{MISC}{
        \itemize{
                \item added vignette folder
        }
    }
}
\section{CHANGES in rsae VERSION 0.1-4 (2012-01-08)}{
    \subsection{NEW FEATURES}{
        \itemize{
            \item added status bar for parametric bootstrap
            \item added arguments 'dec' and 'decor' in 'fitsaemodel.control'
                to choose between different types of decorrelation of the
                residuals and the type of matrix square roots (Cholesky vs. SVD)
        }
    }
    \subsection{MISC}{
        \itemize{
            \item Modified documentation
            \item implemented better checks for the robustness tuning constant
                in '.fitsaemodel.huberm'
        }
    }
    \subsection{BUG FIXES}{
        \itemize{
            \item fixed the names of coefficients returned by 'makedata'
            \item fixed print method 'fitsaemodel' when different robustness
                tuning are used
        }
    }
}
\section{CHANGES in rsae VERSION 0.1-3 (2011-07-26)}{
    \subsection{NEW FEATURES}{
        \itemize{
            \item added function 'convergence' to learn more about the
                convergence of the methods
            \item added argument 'areameans' in function 'robpredict'
            \item added argument 'type' in plot method for class
                'meanssaemodel'
        }
    }
    \subsection{MISC}{
        \itemize{
            \item removed argument 'full' in the summary method for objects
                of class 'fitmodel'
        }
    }
}
\section{CHANGES in rsae VERSION 0.1-2 (2011-07-23)}{
    \itemize{
        \item Initial release
    }
}
