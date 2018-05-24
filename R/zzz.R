
 
.onAttach <- function(libname, pkgname) {
	            
	dataPkg <- try(XML::xmlParse("http://cran.r-project.org/package=everest", isHTML=T), silent=T)
	currVersion <- as.character(utils::packageVersion("everest"))
	
	everest.logo <- paste0("
                                 __         Everest R package 
  ___ _   _____  ________  _____/ /_        ---------------------
 / _ \\ | / / _ \\/ ___/ _ \\/ ___/ __/        LC-MS Metabolite Annotation
/  __/ |/ /  __/ /  /  __/___ / /_
\\___/|___/\\___/_/   \\___/____/\\__/          Version ", currVersion, "

 - Type 'citation('everest')' for citing this R package in publications. 
 - Type 'vignette('everestManual', package='everest')' for a tutorial on Everest's usage. 
 - For bugs, problems and issues, please do not hesitate in contacting xdomingo@scripps.edu.    
    ")
                                      
    packageStartupMessage(everest.logo)   
	
	if(class(dataPkg)[1]!="try-error")
	{
		xml_data <- XML::xmlToList(dataPkg)
		cranVersion <- xml_data$body$table[[1]][[3]]
		#vers_state <- (cranVersion!=currVersion)
		vers_state <- utils::compareVersion(cranVersion, currVersion)
		if(vers_state==-1) vers_state <- 0
		wMsg <- paste("The current version of everest (", currVersion, ") is outdated. There is a new version (", cranVersion, ") available at CRAN. To update your version execute the following: \n install.packages('everest')", sep="")
		if(vers_state) warning(wMsg)
	}else{
		wMsg <- paste("The current available version of everest in CRAN could not be checked. Please, make sure that your current version of everest (", currVersion, ") is the same as in CRAN (http://cran.r-project.org/package=everest). To update your version execute the following: \n install.packages('everest')", sep="")
		warning(wMsg)
	}       	        
	 	
}    



