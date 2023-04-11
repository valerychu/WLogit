\name{CalculPx}
\alias{CalculPx}

\title{
Calculate the class-conditional probabilities.
}
\description{
Calculate the probability for a repsonse to be 1 in the logistic regression model.
}
\usage{
CalculPx(X, beta, intercept = 0)
}
\arguments{
  \item{X}{
Design matrix of the logistic model considered.
}
  \item{beta}{
Vector of coefficients of the logistic model considered.
}
  \item{intercept}{
Whether there is the intercept
}
}

\value{
\item{prob }{the probability for a repsonse to be 1}
}

\seealso{ Please read \url{https://hastie.su.domains/Papers/glmnet.pdf} for more details
}

\author{
Wencan Zhu
}

\examples{
data(X)
data(beta)
CalculPx(X=X, beta=beta)


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, beta, intercept = 0) 
{
    prob <- 1/(1 + exp(-(X \%*\% beta + intercept)))
    return(prob)
  }
}