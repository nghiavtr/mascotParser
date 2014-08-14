\name{MascotParser-class}
\Rdversion{1.1}
\docType{class}
\alias{MascotParser-class}
\alias{getAssignedPeaks,MascotParser-method}
\alias{getKeywordsInfo,MascotParser-method}
\alias{getMSMSqueryInfo,MascotParser-method}
\alias{getPosHits,MascotParser-method}
\alias{getUnassignedPeaks,MascotParser-method}
\alias{parse,MascotParser-method}
\alias{plotSpec,MascotParser-method}
\alias{plotSpecPanel,MascotParser-method}
\alias{setInputfilename,MascotParser-method}

\title{Class \code{"MascotParser"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MascotParser", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Object of class \code{"character"} ~~ }
    \item{\code{remarks}:}{Object of class \code{"character"} ~~ }
    \item{\code{boundaryCode}:}{Object of class \code{"character"} ~~ }
    \item{\code{sourcefile}:}{Object of class \code{"character"} ~~ }
    \item{\code{searchType}:}{Object of class \code{"character"} ~~ }
    \item{\code{fastafile}:}{Object of class \code{"character"} ~~ }
    \item{\code{database}:}{Object of class \code{"character"} ~~ }
    \item{\code{enzyme}:}{Object of class \code{"character"} ~~ }
    \item{\code{taxonomy}:}{Object of class \code{"character"} ~~ }
    \item{\code{header}:}{Object of class \code{"character"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"character"} ~~ }
    \item{\code{qmass}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qexp}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qcharge}:}{Object of class \code{"character"} ~~ }
    \item{\code{qintensity}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qmatch}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qplughole}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{PIHaccessions}:}{Object of class \code{"character"} ~~ }
    \item{\code{PIHhitScores}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{PIHproteinMasses}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{PIHnumHits}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{PIHpeakInfo}:}{Object of class \code{"list"} ~~ }
    \item{\code{MSMSqueryHits}:}{Object of class \code{"list"} ~~ }
    \item{\code{MSMSproteinInfo}:}{Object of class \code{"character"} ~~ }
    \item{\code{MSMSqueryInfo}:}{Object of class \code{"list"} ~~ }
    \item{\code{MSMSnumQuery}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getAssignedPeaks}{\code{signature(object = "MascotParser")}: ... }
    \item{getKeywordsInfo}{\code{signature(object = "MascotParser")}: ... }
    \item{getMSMSqueryInfo}{\code{signature(object = "MascotParser")}: ... }
    \item{getPosHits}{\code{signature(object = "MascotParser")}: ... }
    \item{getUnassignedPeaks}{\code{signature(object = "MascotParser")}: ... }
    \item{parse}{\code{signature(object = "MascotParser")}: ... }
    \item{plotSpec}{\code{signature(object = "MascotParser")}: ... }
    \item{plotSpecPanel}{\code{signature(mz1 = "MascotParser")}: ... }
    \item{setInputfilename}{\code{signature(object = "MascotParser")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MascotParser")
}
\keyword{classes}
