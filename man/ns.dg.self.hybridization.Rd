\name{ns.dg.self.hybridization}
\alias{ns.dg.self.hybridization}
\title{ dG self-hybridization }
\description{
  Calculate the self-hybridization dG for a vector of sequences at a given 
  temerature, using hybrid-ss-min from Unafold package.
}
\usage{
ns.dg.self.hybridization(s, t = 75, c.Na = 1, c.Mg = 0,acid.type=c("DNA","RNA"))
}
\arguments{
  \item{s}{ a vector of sequences in standard IUB/IUPAC nucleic acid codes}
  \item{t}{ the temperature at which dG is calculated }
  \item{c.Na}{ Natrium concentration at which dG is calculated}
  \item{c.Mg}{ Magnesium concentration at which dG is calculated}
  \item{acid.type}{ Calculate dG for DNA or RNA}
}
\details{
  Currently only sequences containing determined nucleotides (ACGT) will produce output.
}
\value{
  A vector of reals giving the dG for all entered sequences
}
\author{ Wim de Leeuw }
\note{ 
  The calculation is performed by calling the external program hybrid-ss-min. This program is part of 
  the UNAFold nucleic acid folding package by Michael Zuker and Nick Markham. This package must 
  be installed to in order to use the function. 
}
\examples{
  ns.dg.self.hybridization(c("ATATATATA","CCCCGGGCGCG","CGCGGGGCATTATATCTACTACTAGTATCT"));
}

