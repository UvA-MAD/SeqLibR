`ns.read.table` <- 
function(ffile,name.col = 1,seq.col=2,header=FALSE)
{
        seq <- read.table(ffile,stringsAsFactors=FALSE,header=header,sep="\t",quote="",comment="")
        seqv <- seq[,seq.col];
        if (name.col > 0)
        {
                names(seqv) <- seq[,name.col]
        }
        seqv;
}

