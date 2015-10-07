\name{ns.read.fasta}
\alias{ns.read.fasta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read fasta file}
\description{
  Read a fasta file and put the read sequences in a vector.
}
\usage{
ns.read.fasta(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ A character vector containing filename to read }
}
\details{
  The function reads a fasta file.
}
\value{
  The function returns a character vector containing the sequences.
  The names of the elements are the set to the comments preceding the sequence.
}
\references{--}
\author{ Wim de Leeuw (w.c.deLeeuw@uva.nl) }
\note{ Charcters not recognized as IUB/IUPAC nucleic acid codes whill trigger a warning. }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }

