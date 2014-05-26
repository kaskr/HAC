\name{[.HAC}
\alias{[.HAC}
\title{Extract tuples.}
\usage{
  [.HAC(x, i, ...)
}
\arguments{
  \item{x}{HAC object}

  \item{i}{Integer vector}
}
\value{
  HAC object
}
\description{
  Extract tuples of HAC object.
}
\details{
  Extract subset of tuples. For instance x[1:2] extracts
  the first two tuples. Alternatively the method can be
  indirectly invoked by the \code{subset} function.
}
\examples{
\dontshow{
hacfile <- system.file("hac","Hac-test_000001.hac",package="HAC")
x <- readExchange(hacfile)
}
x[1:2]
subset(x,type==10000)
split(x,x$type)
}

