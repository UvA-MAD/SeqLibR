\name{ns.blat}
\alias{ns.blat}
\title{ Perform a sequence-alignment using blat}
\description{
  Calculate alignments of one or more query sequences to a 
  sequence database using the blat tool. 
}
\usage{
ns.blat(query, database, eval = 10, tileSize = 10, add.parm = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{ A character vector containing the sequence(s) to be aligned. The character vector is interprented fasta file where elements in the vector are lines. }
  \item{database}{ A reference to the database to align with.}
  \item{eval}{ Numerical threshold e-value. Only alignments with an e-value equal or smaller than \code{eval} will be returned  }
  \item{tileSize}{ tile size which blat will use See blat documentation for details }  
  \item{add.parm}{Additional options to be passed to blat can be set here. }
}
\details{
  The function executes "blat  ....." and retrieves the results. Execute \code{system("blat")} to get an overview of the options which can be passed to the \code{add.parm}. 
}
\value{
  The function returns a data.frame containing the blat result for all submitted sequences; one row per alignment.
  The columns in the output contain the following information with respect to the alignment.
  \item{query}{Name of query sequence.}
  \item{sequence}{Name of sequence in the data base.}
  \item{mlength}{Length of the match.}
  \item{identity}{percent identy of the match.}
  \item{mismatches}{The number of mismatches in the alignment.}
  \item{gaps}{The number of gaps in the alignment.}
  \item{qstart}{The start of the match in the query.}
  \item{qstop}{The end of the match in the query.}
  \item{sstart}{The start of the match in the query.}
  \item{sstop}{The end of the match in the query.}
  \item{eval}{The e-value of the match.}
  \item{bitscore}{The the bit-score of the match.}
}
\author{ Wim de Leeuw (w.c.deLeeuw@uva.nl) }
\examples{
  ts = c("sequence_1"="ATGCGCGTACATCGCCCCCCCGGGGGG","sequence_2"="TCCCCCCCCGGGGGGATCTTATATATATCCCCCGGGGG")
  # Self-self blat 
  ns.blat(ts,ts);
  # Blat some sequecnes
  ns.blat(c("query_A"="TCGCCCCCCCGGGGGG","query_B"="GATCTTATATATATCCC"),ts,my.db,eval=1)
}
