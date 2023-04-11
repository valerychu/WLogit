\name{top_thresh}
\alias{top_thresh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Thresholding to a given threshold of the smallest values}
\description{
This function keeps only the K largest values of the vector \code{vect} and sets the others to the smallest value among the K largest.
}
\usage{
top_thresh(vect, thresh)
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
sorted_vect=sort(x,decreasing=TRUE)
thresh=3
top_thresh(x,thresh)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vect, thresh) 
{
    sorted_vect <- sort(vect, decreasing = TRUE)
    v = sorted_vect[thresh]
    ifelse(vect >= v, vect, v)
  }
}
\keyword{ ~thresholding }