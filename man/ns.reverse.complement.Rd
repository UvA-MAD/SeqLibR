\name{ns.reverse.complement}
\alias{ns.reverse.complement}
\title{ Reverse complement of a vector of sequences}
\description{
  Calculate reverse complement of a sequence
}
\usage{
ns.reverse.complement(s)
}
\arguments{
  \item{s}{ a vector of sequences in standard IUB/IUPAC nucleic acid codes }
}
\details{
  Codes which represent non determined nucleotides will be handled. 
}
\value{
  A vector of sequences giving the reverse complement of the entered sequences
}
\author{ Wim de Leeuw }
\examples{
  ns.reverse.complement(c("ATATATATA","CCCCGGGCGCG","CGCGGGGCATTATATCTACTACTAGTATCT"));
}
