\name{NEWS}
\title{PMCMRplus News}
\encoding{UTF-8}
\section{CHANGES IN PMCMRPLUS VERSION 1.3.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{snkTest} for performing
      the Student-Newman-Keuls test was added.
      \item The funtion \code{duncanTest} for performing
      the Duncan's (new) multiple range test was added.
      \item The funtion \code{welchManyOneTTest} for performing
      many-to-one comparisons with Welch's t-test was added.
      Likewise, \code{welchManyOneTTest} was added as \code{test}
      in function \code{powerMCtests}
      \item The funtion \code{steelTest} including
      print and summary S3-methods for applying
      Steel's Many-to-One Rank Test was added.
      \item The function \code{tamhaneT2Test} was extended
      for performing Tamhane's extended T2' test as well. Thanks to
      Sirio Bolaños.
    }
  }
  \subsection{BUGFIXES}{
    \itemize{
      \item Several Rd files were corrected (e.g. links to other man
      files, corrected calls to functions in examples section)
      \item The function \code{tamhaneT2Test} was corrected in order
      to compute the two-tailed p-value with Sidaks p-adjustment using
      \eqn{m = k(k-1)/2}.
    }
  }
}
\section{CHANGES IN PMCMRPLUS VERSION 1.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{shirleyWilliamsTest} was extended
      for alternatives less and greater using a look-up table
      for critical Williams-t values including
      print and summary methods.
    }
  }
  \subsection{BUGFIXES}{
    \itemize{
      \item Links in Rd file of bwsManyOneTest were corrected.
    }
  }
}
\section{CHANGES IN PMCMRPLUS VERSION 1.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{williamsTest}
      including print and summary methods were added.
    }
  }
}
\section{CHANGES IN PMCMRPLUS VERSION 1.0.1}{
  \subsection{BUGFIXES}{
    \itemize{
      \item Misspelled author name Terpstra was
      corrected in function \code{jonckheereTest}
      and corresponding Rd file.
      \item Reference in the title of Rd files
      for functions \code{mandelhTest} and \code{mandelkTest}
      were corrected.
    }
  }
}
