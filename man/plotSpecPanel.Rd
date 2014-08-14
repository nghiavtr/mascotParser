\name{plotSpecPanel}
\alias{plotSpecPanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add later
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSpecPanel(mz1, int1, mz2, int2, plotTitle = NULL, xlimVal = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mz1}{
%%     ~~Describe \code{mz1} here~~
}
  \item{int1}{
%%     ~~Describe \code{int1} here~~
}
  \item{mz2}{
%%     ~~Describe \code{mz2} here~~
}
  \item{int2}{
%%     ~~Describe \code{int2} here~~
}
  \item{plotTitle}{
%%     ~~Describe \code{plotTitle} here~~
}
  \item{xlimVal}{
%%     ~~Describe \code{xlimVal} here~~
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
structure(function (mz1, int1, mz2, int2, plotTitle = NULL, xlimVal = NULL) 
{
    standardGeneric("plotSpecPanel")
  }, generic = structure("plotSpecPanel", package = "MascotParser"), package = "MascotParser", group = list(), valueClass = character(0), signature = c("mz1", 
"int1", "mz2", "int2", "plotTitle", "xlimVal"), default = `\001NULL\001`, skeleton = (function (mz1, 
    int1, mz2, int2, plotTitle = NULL, xlimVal = NULL) 
stop("invalid call in method dispatch to 'plotSpecPanel' (no default method)", 
    domain = NA))(mz1, int1, mz2, int2, plotTitle, xlimVal), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
