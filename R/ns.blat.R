`ns.blat` <-
function(query,database,eval=10,tileSize=10,add.parm="")
{
	if (!is.character(database)  || (length(database) != 1))
	{
                cat("First initialize with ns.blast.db\n");
                return(FALSE);
	}
	if (!is.null(names(query)) && !ns.blast.unique.names(names(query)))
	{
		cat("sequence names do not have unique first word.");
		return(FALSE)
	}

	eblast.cmd = "blat"
        add.parm = paste(add.parm," -tileSize=",tileSize,sep="");
        seq.file <- tempfile(pattern = "bsequence")
        db.file <- tempfile(pattern = "bdb")
        ns.write.fasta(seq.file,query);
        ns.write.fasta(db.file,database);
        br <- system(paste(eblast.cmd,db.file,seq.file,add.parm,"-out=blast8 stdout"),intern=TRUE)
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
		bt <- bt[bt$eval < eval,]
        } else {
	    bt <- NULL;
	}
	unlink(seq.file);
	bt
}

