\name{getPosHits}
\alias{getPosHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add later
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPosHits(object, scoreThres = 41)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{scoreThres}{
%%     ~~Describe \code{scoreThres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, scoreThres = 41) 
{
    standardGeneric("getPosHits")
  }, generic = structure("getPosHits", package = "MascotParser"), package = "MascotParser", group = list(), valueClass = character(0), signature = c("object", 
"scoreThres"), default = `\001NULL\001`, skeleton = (function (object, 
    scoreThres = 41) 
stop("invalid call in method dispatch to 'getPosHits' (no default method)", 
    domain = NA))(object, scoreThres), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
