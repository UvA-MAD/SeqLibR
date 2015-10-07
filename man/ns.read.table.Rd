\name{ns.read.table}
\alias{ns.read.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read sequence from a tab separated file}
\description{
  Read a tab separated file and put the read sequences in a character vector.
}
\usage{
  ns.read.table(ffile,name.col = 1,seq.col=2,header=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ffile}{ A character vector containing filename to read } 
  \item{seq.col}{ The column number of the table containing the sequence} 
  \item{name.col}{ The column number of the table containing an identifier of the sequence. If the value 0 is given. 
  Then no identifiers are assigned } 
  \item{header}{ Flag if the first line of the file contains a header. The content of the header is ignored. }
}
\details{
  The function reads a tab separated file in which one of the columns 
}
\value{
  The function returns a character vector containing the sequences. 
  The names of the sequences are the set if the name.col to the column number containing the names.
}
\references{--}
\author{ Wim de Leeuw (w.c.deLeeuw@uva.nl) }
\note{ No error checking is done currently so beware.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
