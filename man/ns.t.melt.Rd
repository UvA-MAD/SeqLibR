\name{ns.t.melt}
\alias{ns.t.melt}
\title{ Calculate melting temerature of sequence}
\description{
  Use Nearest neighbor method with parameters by Santa Lucia to estimate melting temperature of a sequence. 
}
\usage{
ns.t.melt(s, ct = 1e-12)
}
\arguments{
  \item{s}{ a vector of sequences in standard IUB/IUPAC nucleic acid codes }
  \item{ct}{ sequence concentration  }
}
\details{
    Sequences with non determined nucleotides (not A,C,G or T) will return NA 
}
\value{
  A vector of reals giving the melting temperature of all entered sequences
}
\author{ Wim de Leeuw }
\examples{
  ns.t.melt(c("ATATATATA","CCCCGGGCGCG","CGCGGGGCATTATATCTACTACTAGTATCT"));
}
