\name{ns.blast.db}
\alias{ns.blast.db}
\title{ Setup a database for blasting sequences }
\description{
Prepare a database for blasting sequences. It either can be an existing database 
or a temporary database for the current R-session only. 
}
\usage{
ns.blast.db(fasta.files = "", sequences = NULL, db.dir = "", blast.path = "")
}
\arguments{
  \item{fasta.files}{An array of files which contain the fasta sequences to be put in the database }
  \item{sequences}{A vector of sequences to put in the database}
  \item{db.dir}{Location where the data base will be stored }
  \item{blast.path}{Location of the blast executable. Use this to parameter to specify the location of the blast tools if the location is not in the PATH-variable}
}
\value{
If succesfull the function will return an database identifier which
can be passed by the \code{database} parameter to the function \code{ns.blast} 
to identify the sequence database to blast against. FALSE will be returned 
if database initialization failed somehow.
}

\details{
either \code{db.dir}, \code{fasta.files} or \code{sequences} has to
be set. If only db.dir is set the previously created database in that
directory is used. This database should be created using the eblast
package. If \code{db.dir} is not set a temporary database is created.
The list of files in fasta.files and the sequences in \code{sequences} are
added to the database. Temporary databases will be removed at the end of
the R session. If \code{db.dir} is set and also  \code{fasta.files} or 
\code{sequences} the given sequences will be added permanently to the 
existing database. 

Each time the function is called without a \code{db.dir} parameter a new temporary 
database is created in the \code{tempdir()}. Reuse the database if blasting many 
times against the same database if you don't want to fill up the disk quickly (See examples).
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
  # better would be ...
  my.db <- ns.blast.db(sequences=ts)
  ns.blast(ts,my.db)
  # ...because we can reuse the database to blast some other sequences
  ns.blast(c("query_A"="TCGCCCCCCCGGGGGG","query_B"="GATCTTATATATATCCC"),my.db,eval=1)
}
