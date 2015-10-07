\name{ns.cg.content}
\alias{ns.cg.content}
\title{ CG-content }
\description{
  Calculate CG content of a sequence
}
\usage{
ns.cg.content(s)
}
\arguments{
  \item{s}{ a vector of sequences in standard IUB/IUPAC nucleic acid codes }
}
\details{
  Codes which represent non determined nucleotides will be countet as 
  if the probability for each nucleotide is equal. For example D (G, A or T) 
  will be counted as 0.33333. 
}
\value{
  A vector of reals giving the cg-content of all entered sequences
}
\author{ Wim de Leeuw }
\examples{
  ns.cg.content(c("ATATATATA","CCCCGGGCGCG","CGCGGGGCATTATATCTACTACTAGTATCT"));
}
