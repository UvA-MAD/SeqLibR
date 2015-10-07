\name{ns.clean}
\alias{ns.clean}
\title{Clean sequences} 
\description{
   Remove linker and adapter sequences from a vector of sequences. 
}
\usage{
ns.clean(seqs,adapters,min.len=60,do.rc=TRUE,eval=1e-5)
}
\arguments{
  \item{seqs}{a vector of sequences in standard IUB/IUPAC nucleic acid codes }
  \item{adapters}{vector containing sequences to cut from seqs} 
  \item{min.len}{Minimum length of output sequences}
  \item{do.rc}{Look for reverse complement of adapters as well} 
  \item{eval}{expectation value threshold of a match to adapter sequence} 
}
\details{
 The function aligns the input sequences to the adapter seqences using Blast and 
removes all matches from the sequences. All remaining parts of the sequences 
which have a length of min.len or more are output. 
}
\value{
A new vector of sequences. The names of sequences are equal the names 
of the input sequences if the sequence is unchanged. If the sequence contains 
the adapter, the output sequence(s) are have the name of the original sequence 
with a part number attached (e.g. EST_ABC with an adapter somewhere in the 
middle will become EST_ABC_p1, EST_ABC_p2).
}
\author{ Wim de Leeuw }
