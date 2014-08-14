\name{QueryHits-class}
\Rdversion{1.1}
\docType{class}
\alias{QueryHits-class}
\alias{addQueryHitInfo,QueryHits-method}
\alias{parseQueryhit,QueryHits-method}

\title{Class \code{"QueryHits"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QueryHits", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{queryHitText}:}{Object of class \code{"character"} ~~ }
    \item{\code{peptideSequence}:}{Object of class \code{"character"} ~~ }
    \item{\code{peptideMass}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{hitScore}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{addQueryHitInfo}{\code{signature(object = "QueryHits")}: ... }
    \item{parseQueryhit}{\code{signature(object = "QueryHits")}: ... }
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
showClass("QueryHits")
}
\keyword{classes}
