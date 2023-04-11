\name{CalculWeight}
\alias{CalculWeight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the weight
}
\description{
Calculate the weight in the penalized weighted-
least-squares problem
}
\usage{
CalculWeight(Px)
}

\arguments{
  \item{Px}{
The vector of estimated probability for each response to be 1.
}
}

\value{
Returns a vector containing the weight.
}


\seealso{ Please read \url{https://hastie.su.domains/Papers/glmnet.pdf} for more details
}

\author{
Wencan Zhu
}

\examples{
data(X)
data(beta)
px <- CalculPx(X=X, beta=beta)
CalculWeight(px)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Px) 
{
    return(Px * (1 - Px))
  }
}