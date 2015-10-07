`ns.read.fasta` <-
function(file) { 
  s <- .Call("seqlib_read_fasta",file,PACKAGE="SeqLibR") 
  if (length(s) != length(unique(names(s)))) 
	warning("Sequence names are not unique.");
  s
}


