`ns.gen.tile` <- 
function(seq,len=60,step=1,circular=0) { 
  .Call("seqlib_gen_tile", seq, as.integer(len), as.integer(step),as.integer(circular),PACKAGE="SeqLibR") 
}

