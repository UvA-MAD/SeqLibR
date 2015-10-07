\name{ns.blast.unique.names}
\alias{ns.blast.unique.names}
\title{ Check sequence names for uniqueness }
\description{
Helper function to see if sequence names can be used in blast alignment.
}
\usage{
ns.blast.unique.names(names)
}
\arguments{
  \item{names}{ vector of sequence names}
}
\details{
Blast returns the first word of the fasta comment as an id. If multiple sequences have the 
same first word it would be impossible to reliable interpret the result. This function can 
be used to check uniqueness.
}
\value{
  TRUE is the names are unique, FALSE otherwise.
}
\author{ Wim de Leeuw }
\examples{
  ns.blast.unique.names(c("sequence1","sequence2","sequence3"))
  ns.blast.unique.names(c("myset sequence1","myset sequence2","myset sequence3"))
	
}
