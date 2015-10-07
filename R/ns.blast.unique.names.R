`ns.blast.unique.names` <- 
function(names) 
{
   length(names)==length(unique(unlist(try(lapply(strsplit(names," "),function(x) { x[[1]] }),silent=TRUE))))
}

