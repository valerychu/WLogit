\name{Thresholding}
\alias{Thresholding}

\title{
Thresholding on a vector
}
\description{
This function provides the thresholding (correction) given a vector. It calls the function \code{\link{top}} or \code{\link{top_thresh}} in the same package, and the output is the vector after correction with the optimal threshold parameter.
}

\usage{
Thresholding(X, y, coef, TOP)
}

\arguments{
  \item{X}{
Design matrix of the logistic model considered.
}
  \item{y}{
Binary response
}
  \item{coef}{
Candidate vector to be corrected
}
  \item{TOP}{
The grill of thresholding
}
}

\value{
\item{opt_top }{The optimal threshold}
\item{auc }{the log-likelihood for each grill of thresholding}
}

\author{
Wencan Zhu
}