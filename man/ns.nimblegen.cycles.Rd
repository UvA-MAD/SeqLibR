\name{ns.nimblegen.cycles}
\alias{ns.nimblegen.cycles}
\title{ Nimblegen Cycles }
\description{
  Calculate the number of cycles which would be required to produce the sequence 
  on a Nimblegen micro-array. 
}
\usage{
ns.nimblegen.cycles(s)
}
\arguments{
  \item{s}{ a vector of strings representing the sequences. }
}
\value{
  A vector of integers giving the number of cycles needed of all entered sequences
}
\author{ Wim de Leeuw }
\examples{
  ns.nimblegen.cycles(c("ATATATATA","CCCCGGGCGCG","CGCGGGGCATTATATCTACTACTAGTATCT"));
}

