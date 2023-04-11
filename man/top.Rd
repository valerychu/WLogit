\name{top}
\alias{top}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Thresholding to zero of the smallest values}
\description{This function keeps only the K largest values of the vector \code{sorted_vect} and sets the others to zero.
}
\usage{
top(vect, thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vect}{
vector to threshold
}
  \item{thresh}{
threshold
}
}
\value{
This function returns the thresholded vector.
}

\author{
Wencan Zhu
}

\examples{
x=sample(1:10,10)
thresh=3
top(x,thresh)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vect, thresh) 
{
    sorted_vect <- sort(abs(vect), decreasing = TRUE)
    v = sorted_vect[thresh]
    ifelse(abs(vect) >= v, vect, 0)
  }
}
\keyword{ ~thresholding }
