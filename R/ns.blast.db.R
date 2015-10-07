`ns.blast.db` <-
function (fasta.files="",sequences=NULL,db.dir="",blast.path="") {
  check.exe <- function (name) { if (!system(paste("which ",name,"> /dev/null",sep=""),ignore.stderr = TRUE)) TRUE else FALSE }
	blast.cmd = "blastn";
	format.cmd = "makeblastdb";
	if (blast.path!="") {
		blast.cmd =paste(as.character(blast.path),"/",blast.cmd,sep="")
		format.cmd =paste(as.character(blast.path),"/",format.cmd,sep="")
	}
	if (!check.exe(format.cmd) || !check.exe(blast.cmd))  {
		cat("Blast intialization FAILED: Blast execuatables where not found\n")
		return(FALSE)
	}
	if ((db.dir == "")  && (fasta.files=="") && (is.null(sequences))) {
		cat("Blast intialization FAILED: no fasta.files, db.dir or sequences given\n")
		return(FALSE)
	}
	if (db.dir == "") {
		db.dir <- paste(tempdir(),"/BlastDB_",sep="")
		i = 0;
		while (!file.access(paste(db.dir,i,sep=""),mode=0))  {
			i = i +1;
		}
		db.dir = paste(db.dir,i,sep="")
		system(paste("mkdir",db.dir))
		if (!file.info(db.dir)$isdir) {
			cat("Blast intialization FAILED: could not create db dir\n")
			return(FALSE)
		}
	}
	if (is.na(file.info(db.dir)$isdir)) {
		cat(paste("Creating blast database directory : ",db.dir,"\n"));
		system(paste("mkdir",db.dir))
	}
	if (is.na(file.info(db.dir)$isdir) | !file.info(db.dir)$isdir) {
		cat(paste("Blast intialization FAILED could not access db.dir :",db.dir,"\n"))
		return(FALSE)
	}
	if (fasta.files[1] !="") {
		for (ff in as.character(fasta.files)) {
			if (file.access(ff,mode=4) == 0) {
				if (system(paste("cp ",ff," ",db.dir,"/",basename(ff),".nt; cd ",db.dir,";",format.cmd," -logfile /dev/null -dbtype nucl -in ",basename(ff),".nt",sep="")))
				{
					cat(paste("Blast intialization FAILED could not format file:",ff," as blast database\n"))
					return(FALSE)
				}
			}
			else {
				cat(paste("Blast intialization FAILED could not find sequence file:",ff,"\n"))
				return(FALSE)
			}

		}
	}
	if (!is.null(sequences))
	{
		if (!is.character(sequences))
		{
			cat("sequences variable does not contain charater data");
			return(FALSE)
		}
		if (!is.null(names(sequences)) && !ns.blast.unique.names(names(sequences)))
		{
			cat("sequence names do not have unique first word.");
			return(FALSE)
		} 
		ff = tempfile("bseqs");
		ns.write.fasta(ff,sequences);
		if (system(paste("cp ",ff," ",db.dir,"/",basename(ff),".nt; cd ",db.dir,";",format.cmd," -logfile /dev/null -dbtype nucl -in ",basename(ff),".nt",sep="")))
		{
			cat(paste("Blast intialization FAILED could not format sequences as blast database\n"))
			return(FALSE)
		}
	}
	dbs <- paste(system(paste("ls -1 ",db.dir,"/*.nt",sep=""),intern=TRUE),sep="",collapse=" ")
	if (dbs == "")
	{
		cat(paste("Blast intialization FAILED the specified directory is empty\n"))
		return(FALSE)
	}
	paste(blast.cmd,' -outfmt 6 -db "',dbs,'" ',sep="");
}

