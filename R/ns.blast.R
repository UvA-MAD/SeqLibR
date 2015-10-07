`ns.blast` <-
function(x,database,eval=10,filter=FALSE,add.parm="")
{
	if (!is.character(database)  || (length(database) != 1))
	{
                cat("First initialize with ns.blast.db\n");
                return(FALSE);
	}
	if (!is.null(names(x)) && !ns.blast.unique.names(names(x)))
	{
		cat("sequence names do not have unique first word.");
		return(FALSE)
	}

	eblast.cmd = database;
        if (!filter) add.parm=paste(add.parm,"-dust no");
        add.parm = paste(add.parm,"-evalue",eval);
        seq.file <- tempfile(pattern = "bsequence")
        ns.write.fasta(seq.file,x);
        br <- system(paste(eblast.cmd,add.parm,"-query",seq.file),intern=TRUE)
	k <- unlist(strsplit(br,"\t"))
	klen = length(k)
	if (klen > 0) {
	    bt <- data.frame(
		query    = as.character(k[seq(1,klen,by=12)]),
		sequence = as.character(k[seq(2,klen,by=12)]),
		identity =   as.numeric(k[seq(3,klen,by=12)]),
		mlength  =   as.integer(k[seq(4,klen,by=12)]),
		mismatches=  as.integer(k[seq(5,klen,by=12)]),
		gaps      =  as.integer(k[seq(6,klen,by=12)]),
		qstart    =  as.integer(k[seq(7,klen,by=12)]),
		qstop     =  as.integer(k[seq(8,klen,by=12)]),
		sstart    =  as.integer(k[seq(9,klen,by=12)]),
		sstop     =  as.integer(k[seq(10,klen,by=12)]),
		eval      =  as.numeric(k[seq(11,klen,by=12)]),
		bitscore  =  as.numeric(k[seq(12,klen,by=12)]),stringsAsFactors=FALSE)
        } else {
	    bt <- NULL;
	}
	unlink(seq.file);
	bt
}

