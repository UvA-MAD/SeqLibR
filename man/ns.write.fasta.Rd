\name{ns.write.fasta}
\alias{ns.write.fasta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write fasta file}
\description{
  Write a fasta file and put the read sequences in a vector.
}
\usage{
ns.write.fasta(ffile,seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ffile}{ The name of the file to write} 
  \item{seqs}{A character vector containing the sequences to write. The element names will be written as fasta comments.} 
}
\details{
 If the elements have names these will be written as fasta comments. If the vector does not contain names the fasta comment will be set to sequence\_\%d, where \%d is the array index 
}
\value{
  The function returns the result of the write.
}
\references{--}
\author{ Wim de Leeuw (w.c.deLeeuw@uva.nl) }
\note{ No error checking is done currently so basically anything will be written but the result might be garbage. }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
