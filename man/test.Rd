\name{test}
\alias{test}
\docType{data}
\title{
WLogit output
}
\description{
The output of WLogit in the example given in the vignette. 
}
\usage{data("test")}
\format{
  The format is:
List of 4
 $ beta          : num [1:50, 1:500] 0 0 0 0 0 ...
 $ lambda        : num [1:50] 100.8 80 73 58.9 56.7 ...
 $ beta.min      : num [1:500] 0.0194 0.0348 0.0259 0.0287 0.0385 ...
 $ log.likelihood: num [1:50] 57.7 57.7 57.7 57.7 57.7 ...
}

\examples{
data(test)
str(test)
}
\keyword{datasets}
