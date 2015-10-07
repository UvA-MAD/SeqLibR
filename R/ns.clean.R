`ns.clean` <-
function(seqs,adapters,min.len=60,do.rc=TRUE,eval=1e-5) { 
	if (is.null(names(seqs)))
	{
		names(seqs) <- paste("sequence_",1:length(seqs),sep="")
	}
	if (!ns.blast.unique.names(names(seqs)))
	{
		warning("sequence names do not have unique first word. Can not blast.");
		return(FALSE)
	}
	dbi = ns.blast.db(sequences=adapters)
	blast.parms=""
	if (!do.rc) { blast.parms = "-S 1" }
	br <- ns.blast(seqs,dbi,eval,filter=FALSE,add.parm=blast.parms);
	if (!is.null(br))
	{
		br$qt <- as.numeric(factor(br$query, levels = factor(unlist(lapply(strsplit(names(seqs), " "), function(x) { x[[1]] })))))
		.Call("sequence_clean",seqs,as.integer(br$qt),as.integer(br$qstart),as.integer(br$qstop),as.integer(min.len),PACKAGE="SeqLibR")
	} else {
		seqs;
	}
}
