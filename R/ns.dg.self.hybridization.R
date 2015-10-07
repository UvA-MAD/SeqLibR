`ns.dg.self.hybridization` <-
function(s,t=75,c.Na=1,c.Mg=0,acid.type= c("DNA","RNA")) {
  tp = match.arg(acid.type);
 .Call("seqlib_dg_selfhyb",s,t,c.Na,c.Mg,as.integer(ifelse(tp=="DNA",1,0)),PACKAGE="SeqLibR") 
}
