\name{ns.gen.tile}
\alias{ns.gen.tile}
\title{ Generate sequence tiles}
\description{
  Given a sequence generate a vector of subsequences 
}
\usage{
  ns.gen.tile(seq,len=60,step=1,circular=0)
}
\arguments{
  \item{seq}{ A sequences in standard IUB/IUPAC nucleic acid codes }
  \item{len}{ the length of the subsequence}
  \item{step}{ the step size }
  \item{circular}{ if not zero assume a circular sequence }
}
\details{
  It is ment for genetic sequences but will work for any string.
}
\value{
  A vector of strings containing the subsequences
}
\author{ Wim de Leeuw }
\examples{
  ns.gen.tile("CGCGGGGCATTATATCTACTACTAGTATCT",10,2);
}
