\name{INSTALL}
\title{How to install \pkg{ltable}}


\section{Installation of v.2.0.1 (2021-08-21) for UNIX(MacOS)}{
  \itemize{
    \item {FIRST install GSL: GNU Scientific Library from  http://www.gnu.org/software/gsl}
    \item {Do regular installation}
    }
}
\section{If problems arise with finding GSL libraries and headers location:}{
  \itemize{
    \item {download the package archive}
    \item {unpack it}
    \item {in shell (can use Terminal in RStudio) type: \cr
  gsl-config --libs \cr
  to get the path to GSL's lib directory (-L<path-to-lib>) \cr
  gsl-config --libs \cr
  to get the path to GSL's lib directory (-L<path-to-lib>)} 
    \item {Change file Makevars in ltable/src:\cr
    Add -L<path-to-lib> to PKG_LIBS \cr
    Add (new) line: PKG_CPPFLAGS=-I<path-to-include>}
    \item {Install ltable via \cr
    LDFLAGS=-L<path-to-lib>; export LDFLAGS \cr
    CPPFLAGS=-I<path-to-include>;export CPPFLAGS
}
    }
}

\section{Installing ltable package standard way may not preserve vignette. To secure vignette :}{
  \itemize{
    \item {download the package archive}
    \item {unpack it, find ltable.pdf in ltable\\vignettes and ltable\\inst\\Doc}
    \item {use RStudio to install: \cr
        - in top right corner choose Project: \cr
        - in pop-up manue choose Open Project...\cr
        - choose package directory \cr
        - in top right menue choose Build \cr
        - in Build submenu choose Install and Restart
    }
}



