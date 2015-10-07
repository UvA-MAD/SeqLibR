\name{ns.blast}
\alias{ns.blast}
\title{ Perform a sequence-alignment }
\description{
  Calculate alignments of one or more sequences in fasta format to a 
  sequence database using the blast tools. The sequence database 
  is set up using init.blast. 
}
\usage{
ns.blast(x, database, eval = 10, filter = FALSE, add.parm = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A character vector containing the sequence(s) to be blasted. The character vector is interprented a fasta file where elements in the vector are lines. }
  \item{database}{ A reference to the database to blast against, use the output of the function \code{ns.blast.db} here } 
  \item{eval}{ Numerical threshold e-value. Only alignments with an e-value equal or smaller than \code{eval} will be returned  }
  \item{filter}{ If \code{filter} is set then repetative sequences in the database and sequence are not considered for alignment. See blastn documentation for details }  
  \item{add.parm}{Additional options to be passed to blastn can be set here. For example to lower the gap penalty use add.parm = "-gapextend -2" }
}
\details{
  The function executes "blastn  ....." and retrieves the results. Execute \code{system("blastn -help")} to get an overview of the options which can be passed to the \code{add.parm}. 
}
\value{
  The function returns a data.frame containing the blast result for all submitted sequences; one row per alignment.
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
\references{Blast: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaeffer,
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),
"Gapped BLAST and PSI-BLAST: a new generation of protein database search
programs",  Nucleic Acids Res. 25:3389-3402.}
\author{ Wim de Leeuw (w.c.deLeeuw@uva.nl) }
\examples{
  ts = c("sequence_1"="ATGCGCGTACATCGCCCCCCCGGGGGG","sequence_2"="TCCCCCCCCGGGGGGATCTTATATATATCCCCCGGGGG")
  # Self-self blast 
  ns.blast(ts,ns.blast.db(sequences=ts))
  # Blast some sequecnes
  my.db <- ns.blast.db(sequences=ts)
  ns.blast(c("query_A"="TCGCCCCCCCGGGGGG","query_B"="GATCTTATATATATCCC"),my.db,eval=1)
}
