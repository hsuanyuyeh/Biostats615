#' @export
recodqual <-
function(X,rename.level=FALSE)
	{
		#X <- as.matrix(X)
		GNA <- tab.disjonctif.NA(X,rename.level)
		G <- replace(GNA,is.na(GNA),0)
		n <- nrow(GNA)
		if (n > 1)
		{
		  ns <- applym(G,2,"sum")
		  nmiss <- applym((isNA(GNA)),2,"sum")
		  if(sum((n-nmiss)==ns)!=0) stop("There are columns in X.quali where all the categories are identical",call.=FALSE)
		}
		return(G)	
	}

