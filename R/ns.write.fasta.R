`ns.write.fasta` <-
function(ffile,seqs)
{
        if (is.null(names(seqs))) {
                write(file=ffile,paste(">sequence_",1:length(seqs),"\n",seqs,sep=""))
        } else {
                write(file=ffile,paste(">",names(seqs),"\n",seqs,sep=""))
        }

}

