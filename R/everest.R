globalVariables(c("hmdb", "everestRules"))

## MetaboSet Class Definition:

setClass(Class = "MetaData", representation = representation(Instrumental = "data.frame", Phenotype = "data.frame", DataDirectory="character"))

setClass(Class = "Statistics", representation = representation(Univariate="data.frame", Multivariate="data.frame"))	

setClass(Class="MSResultsParameters", representation = representation(Alignment = "list", Annotation = "list", Identification = "list"))

setClass(Class="Data", representation = representation(FeatureList = "list", FactorList = "list", Parameters = "list"))

setClass(Class = "Results", representation = representation(Parameters="MSResultsParameters", Alignment = "data.frame", Annotation = "data.frame", Identification="data.frame", Statistics="Statistics"))

setClass(Class="MetaboSet",representation= representation(Info = "character", Data="Data", MetaData="MetaData", Results = "Results"))

setClass(Class = "metaDB", representation = representation(name="character", version="character", info="character", database="list"))

setClass(Class = "evGroup", representation = representation(AlignID="numeric", rt="numeric", features="vector", pseudospectra="character", annotation="data.frame", identification="data.frame"))

setMethod("show", "MetaboSet", function(object){
	Nsamples <- max(c(length(object@Data@FactorList), length(object@Data@FeatureList$data.table)), na.rm=TRUE) 
	cat("A \"MetaboSet\" object containing", length(object@Data@FactorList), "samples \n \n" ,sep=" ")
	cat("Data processed with", object@Data@Parameters$algorithm, "\n" ,sep=" ")
	cat("Info attached to this experiment: \n", object@Info)
})

setMethod("show", "evGroup", function(object){
	cat("An \"Everest Group\" object containing", length(object@features), "features at RT =", object@rt, " \n \n" ,sep=" ")
	nrowMax <- 5
	if(nrowMax>nrow(object@annotation)) nrowMax <- nrow(object@annotation)
	cat("Previsualization of annotation: \n \n")
	print(object@annotation[1:nrowMax,])
	if(nrowMax!=nrow(object@annotation))  cat("... (the rest of the data is not shown)")
	nrowMax <- 5
	if(nrowMax>nrow(object@identification)) nrowMax <- nrow(object@identification)
	cat("\n \n Previsualization of identification: \n \n ")
	if(nrowMax==0) cat('No identification')
	if(nrowMax!=0) {
		print(object@identification[1:nrowMax,])
		if(nrowMax!=nrow(object@identification))  cat("... (the rest of the data is not shown)")
	}
})

metaData <- function(object) {object@MetaData@Instrumental}
phenoData <- function(object) {object@MetaData@Phenotype}
		

paste.sp <- function(x, y) 
{
    paste(x, y, sep = ",")
}

normalize <- function(x) 
{
    x[is.na(x)] <- 0
    if(is.matrix(x) == T) norm.x <- sweep(x, 2, apply(x, 2, function(k) max(k, na.rm = T)), "/")
    if(is.matrix(x) == F) norm.x <- x/max(x, na.rm = T)
    norm.x[is.na(norm.x)] <- 0
    norm.x
}


evAnnotate <- function(xcmsSet=NULL, data.table=NULL, ion.mode=c('pos','neg'), min.correlation=0.7, max.time.dist=1, ppm.error=20, max.mz.error=0.05, annRules=everestRules, db.object=hmdb, fluidType=NULL, topAdductsN=4, import.profiles=FALSE, alternative.data.path=NULL)
{
	if(length(ion.mode)!=1) warning('Ionization mode not selected, taking POSITIVE as default')
	ion.mode <- match.arg(ion.mode, c('pos','neg'))
	annRul <- annRules[which(annRules[,'mode']==ion.mode),]	
	if(nrow(annRul[which(annRul[,'charge']!=1),])!=0)
	{
		warning('The annotation rules selected contain adducts with charge different to 1. Currently, charges different to 1 are not supported. Those have been removed from the annotation rules list.')
		annRul <- annRul[-which(annRul[,'charge']!=1),]	
	}	
		
	if(is.null(xcmsSet) & is.null(data.table)) stop('You must provide with an xcmsSet (xset3) object from XCMS or a data table from other softwares')
	if(!is.null(xcmsSet) & !is.null(data.table)) stop('Please, choose between an xcmsSet object or a data table')

	if(!is.null(xcmsSet))
	{
		## Check if XCMS is installed
		if(!requireNamespace("xcms", quietly = TRUE)) stop("\n \t XCMS is not installed. To use Everest with xcmsSet from XCMS you need to install XCMS first. To install the XCMS package, please visit the bioconductor website: http://bioconductor.org/packages/release/bioc/html/xcms.html\nOr, alternatively, execute the following R code:\n\t\t\n\t\t## try http:// if https:// URLs are not supported \n\t\tsource('https://bioconductor.org/biocLite.R')\n\t\tbiocLite('xcms')")
	
		Experiment <- xcms2MetaboSet(xcmsSet)
		cat('Extracting pseudo-spectra... \n')
		al.list <- groupIonsXCMS(Experiment@Data@FeatureList[[1]], min.lin.rel=min.correlation, max.time.dist=max.time.dist)
		cat('\n A total of', nrow(al.list), 'pseudo-spectra found across samples')
		cat('\n Annotating pseudo-spectra... \n ')
		an.list <- annotFeats(al.list, annRul=annRul, ppm.error=ppm.error, max.mz.error=max.mz.error, db.object=db.object, fluidType=fluidType, topAdductsN=topAdductsN)
		rtMat <- cbind(xcms::groupnames(Experiment@Data@FeatureList[[1]]),  Experiment@Data@FeatureList[[1]]@groups[,c('rtmed')]) 
		Rt <- sapply(an.list[,'Feature'], function(x) rtMat[which(rtMat[,1] %in% x),2])
		anList <- cbind(an.list[,c('AlignID','AnnID', 'Feature', 'Mass')], Rt, an.list[,c('Adduct', 'IsoCount', 'Score', 'mIsoMass', 'Isotope','toMSMS')])
		rownames(anList) <- NULL	
		cat('\n Done!')
	}
	if(!is.null(data.table))
	{	
		if(!any(colnames(data.table)[2:3] == c("mz","RT"))) stop('Your table is not formatted as it should be. Please format your table accordingly. Please see the documentation for more information. Type vignette(',"'",'everestManual',"'",'package=',"'",'erah',"'",')')
		Experiment <- table2MetaboSet(data.table)
		cat('Extracting pseudo-spectra... \n')
		al.list <- groupIonsTable(data.table, min.lin.rel=min.correlation, max.time.dist=max.time.dist)
		cat('\n A total of', nrow(al.list), 'pseudo-spectra found across samples')
		cat('\n Annotating pseudo-spectra... \n ')
		an.list <- annotFeats(al.list, annRul=annRul, ppm.error=ppm.error, max.mz.error=max.mz.error, db.object=db.object, fluidType=fluidType, topAdductsN=topAdductsN)
		rtMat <- cbind(as.character(data.table[,1]),  data.table[,3 ]) 
		Rt <- sapply(an.list[,'Feature'], function(x) rtMat[which(rtMat[,1] == x),2])
		anList <- cbind(an.list[,c('AlignID','AnnID', 'Feature', 'Mass')], Rt, an.list[,c('Adduct', 'IsoCount', 'Score', 'mIsoMass','Isotope','toMSMS')])
		rownames(anList) <- NULL	
		cat('\n Done!')		
	}
	if(is.null(db.object)) db.Name <- NULL
	if(!is.null(db.object)) db.Name <- db.object@name
	Experiment@Results@Parameters@Alignment <- list(algorithm='Everest-ion-grouping', min.correlation=min.correlation, max.time.dist=max.time.dist)
	Experiment@Results@Parameters@Annotation <- list(algorithm='Everest', annRul=annRul, ppm.error=ppm.error, max.mz.error=max.mz.error, db.object.name=db.Name, fluidType=fluidType, topAdductsN=topAdductsN)
	Experiment@Results@Alignment <- al.list	
	Experiment@Results@Annotation <- anList	
	Experiment
}

evIdentify <- function(Experiment, db=hmdb, ppm.error=5){
	#Experiment <- ex
	#db=hmdb
	#ppm.error=10
	
	
	# AlIds <- unique(as.numeric(as.vector(Experiment@Results@Annotation$AlignID[which(as.vector(Experiment@Results@Annotation$Feature)=='M134T107')])))
	# groupAnList <- Experiment@Results@Annotation[which(as.numeric(as.vector(Experiment@Results@Annotation$AlignID))==AlIds),]
	# i <- 8795
	
	misovect <- as.numeric(as.vector(unlist(lapply(db@database, function(x) {misom <- x$monisotopic_molecular_weight
		if(is.null(misom)) misom <- 0
		misom
		}))))
	idTab <- matrix(0, nrow=0, ncol=8)
	cat('Comparing monositopic masses with library...')
	for(i in 1:nrow(Experiment@Results@Annotation))
	{
		qMass <- as.numeric(as.vector(Experiment@Results@Annotation$mIsoMass[i]))
		maxMassDiff <- getMassError(qMass, ppm.error)
		fputCand <- which(abs(qMass - misovect)<maxMassDiff)
		if(length(fputCand)==0) next
		ppmev <- sapply(abs(qMass - misovect[fputCand]), function(x) getPPM(x, qMass))
		
		outT <- sapply(1:length(ppmev), function(x) matrix(c(as.matrix(Experiment@Results@Annotation[i,c('AlignID','AnnID','Feature')]), fputCand[x], db@database[[fputCand[x]]]$accession, db@database[[fputCand[x]]]$name, round(as.numeric(as.vector(db@database[[fputCand[x]]]$monisotopic_molecular_weight)),5), round(ppmev[x],2)), ncol=8), simplify=F)
		outT <- do.call('rbind', outT)
		
		idTab <- rbind(idTab, outT)			
	}
	
	idTabf <- as.data.frame(idTab)
	colnames(idTabf) <- c('AlignID','AnnID','Feature','DB.Id','Accession','Name','mIsoMass(DB)','ppm.error')
		
	Experiment@Results@Parameters@Identification <- list(algorithm='Everest', ppm.error=ppm.error,db=paste(db@name, 'version', db@version))
	Experiment@Results@Identification <- idTabf
	Experiment
}


showAn <- function(Experiment, feature)
{
	AlIds <- unique(as.numeric(as.vector(Experiment@Results@Annotation$AlignID[which(as.vector(Experiment@Results@Annotation$Feature)==feature)])))
	groupAnList <- Experiment@Results@Annotation[which(as.numeric(as.vector(Experiment@Results@Annotation$AlignID))==AlIds),]
	annIds <- as.numeric(as.vector(groupAnList[which(groupAnList[,'Feature'] %in% feature), 'AnnID']))
	groupAnList[which(groupAnList$AnnID %in% annIds),]
}

showGroup <- function(Experiment, feature)
{
	AlIds <- unique(as.numeric(as.vector(Experiment@Results@Alignment$AlignID[grep(feature, as.vector(Experiment@Results@Alignment$Feature))])))
	alTable <- Experiment@Results@Alignment[which(as.numeric(as.vector(Experiment@Results@Alignment$AlignID))==AlIds),]
	anTable <- Experiment@Results@Annotation[which(as.numeric(as.vector(Experiment@Results@Annotation$AlignID))==AlIds),]
	idTable <- Experiment@Results@Identification[which(as.numeric(as.vector(Experiment@Results@Identification$AlignID))==AlIds),]
	
	feat.vct <- strsplit(as.character(alTable$Features), ",")[[1]]

	outObj <- new("evGroup",  AlignID=AlIds, rt=alTable$RT, features=feat.vct, pseudospectra=as.character(alTable$Spectra), annotation=anTable, identification=idTable)
		
	outObj
}

# showId <- function(Experiment, feature)
# {
	# locIdTab <- Experiment@Results@Identification[which(as.vector(Experiment@Results@Identification$Feature)==feature),]
# }

# xmcs2everest <- function(xcmSet, al.path=NULL)
# {
	
	# #xcmSet <- xset3
	# #al.path <- 'Samples'
	
	# ### Pasar les mostres
	
	# if(is.null(al.path)) filePaths <- xcms::filepaths(xcmSet)
	# if(!is.null(al.path)) {
		# al.path.split <- strsplit(al.path, "")[[1]]
		# if(al.path.split[length(al.path.split)] != "/") al.path <- paste(al.path, "/", sep="") 
		# fileNames <- as.vector(unlist(sapply(xcms::filepaths(xcmSet), function(x){
			# y <- strsplit(x, "/")[[1]] 
			# #paste(al.path, y[length(y)], sep="")
			# y[length(y)]
		# })))
		# fileClass <- xcmSet@phenoData
		# filePaths <- as.vector(apply(cbind(fileClass, fileNames),1, function(x) paste(al.path, x[1], '/', x[2], sep='')))
	# } 

	# peakList <- cbind(xcmSet@peaks, rep("",nrow(xcmSet@peaks)))
	# colnames(peakList)[ncol(peakList)] <- "profile"
	# groupList <- xcmSet@groups
	
	# #filepaths(xcmSet@xcmsSet) <- filePaths
	# xcms::filepaths(xcmSet) <- filePaths

	# profileList <- vector('list',length(filePaths))
	# names(profileList) <- rownames(xcmSet@phenoData)
	
	# pb <- txtProgressBar(min=1,max=length(filePaths), width=50, style=3)	
	# for(nSample in 1:length(filePaths))
	# {
		# nFeats <- nrow(xcmSet@groups)
		# tmp <- retEIC(xcmSet, filePaths,index=rep(nSample,nFeats))
		# RTvect <- unlist(tmp$scantimes)
		# profileList[[nSample]] <- vector('list',nFeats)
		# for(nPeak in 1:nFeats)
		# {
			# sigInds <- which(!is.na(tmp$EIC[nPeak,]))
			# profile.text <- paste(sweep(as.matrix(RTvect[sigInds]),1,as.matrix(tmp$EIC[nPeak,sigInds]),"paste.sp"), collapse=" ")
			# profileList[[nSample]][[nPeak]] <- profile.text
		# }
		# setTxtProgressBar(pb, nSample)	
	# }
	# profileList
# }	
		

xcms2MetaboSet <- function(xcmSet)
{	
	
	ident.list <- as.data.frame(matrix(0,ncol=7, dimnames=list(row=0,col= c("AlignID", "tmean", "Name", "MatchFactor", "CAS", "Formula", "DB.Id"))))
	uni.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("AlignID", "FoldChangue", "pvalue"))))
	multi.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("AlignID", "CompoundsInvolved", "pvalue"))))
	
	stat.parameters <- new("MSResultsParameters", Alignment=list(), Identification=list())
	statistics <- new("Statistics", Univariate = uni.stats, Multivariate = multi.stats)
	MS.Results <- new("Results", Parameters = stat.parameters, Identification = ident.list, Statistics = statistics )
		
	### Instrumental /Phenotype conversion:
		
	filePaths <- xcms::filepaths(xcmSet)
	sepList <- strsplit(filePaths, "/")
	sepMat <- do.call("cbind", sepList)
	commonInds <- which(apply(sepMat, 1, function(x) all(x == x[1])))
	path.dir <- paste(sepList[[1]][commonInds], collapse="/")	
	fileNames <- as.vector(unlist(lapply(sepList, function(x) paste(x[-commonInds], collapse="/"))))
	sampleNames <- as.vector(unlist(lapply(sepList, function(x) strsplit(x[length(x)], "\\.")[[1]][1])))
	exClasses <- as.vector(unlist(lapply(sepList, function(x) x[-commonInds][1])))
	
	instrumental.dataframe <- as.data.frame(cbind(sampleID=sampleNames, filename=fileNames))
	phenotype.dataframe <- as.data.frame(cbind(sampleID=sampleNames, class=exClasses))
	
	MS.MetaData <- new("MetaData", Instrumental = instrumental.dataframe, Phenotype = phenotype.dataframe, DataDirectory=path.dir)

	### Peak List to Feature List conversion:
	
	feature.list <- list(xcmsSet=xcmSet)
	#feature.list <- lapply(unique(xcmSet@peaks[,"sample"]), function(x) as.data.frame(xcmSet@peaks[which(xcmSet@peaks[,"sample"]==x),]))		
	#names(feature.list) <- as.vector(instrumental.dataframe$sampleID)
		
	align.list <- list(peakExternal=list(values=matrix(), idgroup=list()), everest=data.frame())
		
	#xcmA <- groupval(xcmSet, value="intb") 
	#xcmI <- groupval(xcmSet, value="maxo") 
	#XCMt <- data.frame(xcmSet@groups)
	#rownames(XCMt)=groupnames(xcmSet)
	#featureTable <- cbind(XCMt$mzmed,XCMt$rtmed, xcmA)	
	#colnames(featureTable)[1:2] <- c("mz","RT")
	
	#align.list$peakExternal$value <-  featureTable
	#align.list$peakExternal$idgroup <- xcmSet@groupidx	
	#align.list$everest <- list()
	
	MS.Data <- new("Data", FeatureList = feature.list, FactorList = list(), Parameters = list(NULL))
	info <- paste(utils::capture.output(xcmSet), collapse='\n')
	sample.container <- new("MetaboSet", Info = info, Data = MS.Data, MetaData = MS.MetaData, Results = MS.Results)
	sample.container@Data@Parameters <- list(algorithm="XCMS", parameters="unknown")
	#sample.container@Results@Parameters@Alignment <- list(algorithm="everest"
	sample.container	
}

table2MetaboSet <- function(data.table)
{	
	
	ident.list <- as.data.frame(matrix(0,ncol=7, dimnames=list(row=0,col= c("AlignID", "tmean", "Name", "MatchFactor", "CAS", "Formula", "DB.Id"))))
	uni.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("AlignID", "FoldChangue", "pvalue"))))
	multi.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("AlignID", "CompoundsInvolved", "pvalue"))))
	
	stat.parameters <- new("MSResultsParameters", Alignment=list(), Identification=list())
	statistics <- new("Statistics", Univariate = uni.stats, Multivariate = multi.stats)
	MS.Results <- new("Results", Parameters = stat.parameters, Identification = ident.list, Statistics = statistics )
		
	### Instrumental /Phenotype conversion:
		
	path.dir <- ''
	fileNames <- rep('Unknown', c(ncol(data.table)-3))
	sampleNames <- as.vector(colnames(data.table[,-c(1:3)]))
	exClasses <- rep('Unk', c(ncol(data.table)-3))
	
	instrumental.dataframe <- as.data.frame(cbind(sampleID=sampleNames, filename=fileNames))
	phenotype.dataframe <- as.data.frame(cbind(sampleID=sampleNames, class=exClasses))
	
	MS.MetaData <- new("MetaData", Instrumental = instrumental.dataframe, Phenotype = phenotype.dataframe, DataDirectory=path.dir)

	### Peak List to Feature List conversion:
	
	feature.list <- list(data.table=data.table)	
	align.list <- list(peakExternal=list(values=matrix(), idgroup=list()), everest=data.frame())
			
	MS.Data <- new("Data", FeatureList = feature.list, FactorList = list(), Parameters = list(NULL))
	info='No info'
	sample.container <- new("MetaboSet", Info = info, Data = MS.Data, MetaData = MS.MetaData, Results = MS.Results)
	sample.container@Data@Parameters <- list(algorithm="Third-party(data table)", parameters="unknown")
	#sample.container@Results@Parameters@Alignment <- list(algorithm="everest"
	sample.container	
}

annoTable <- function(Experiment)
{
	if(Experiment@Data@Parameters$algorithm=='XCMS'){
		fT <- data.frame(Experiment@Data@FeatureList$xcmsSet@groups)
        feat.names <- xcms::groupnames(Experiment@Data@FeatureList$xcmsSet) 
		#rownames(fT) <- feat.names
	}
	if(Experiment@Data@Parameters$algorithm=='Third-party(data table)'){
		fT <- Experiment@Data@FeatureList$data.table
		feat.names <- fT[,1] 
		#rownames(fT) <- feat.names
	}
	
	fT <- cbind(fT, matrix('', nrow=nrow(fT), ncol=5), stringsAsFactors=FALSE)
	colnames(fT)[(ncol(fT) - 4):ncol(fT)] <- c('Adduct', 'AnnID', 'AlignID','Isotope','toMSMS')
	fT$toMSMS <- factor(fT$toMSMS, levels=c(levels(fT$toMSMS), c('yes',''))) 
	fT$Isotope <- factor(fT$Isotope, levels=c(levels(fT$Isotope), c('yes',''))) 
	
	alignIDTab <- apply(Experiment@Results@Alignment, 1, function(x){
		#cat(x, '\n')
		featNames <- strsplit(as.character(x['Features']), ',')[[1]]
		retX <- cbind(featNames, rep(x['AlignID'], length(featNames)))
		retX
	})
	alignIDTab <- do.call('rbind', alignIDTab)
	
	for(i in 1:length(feat.names))
	{
		subAnn <- Experiment@Results@Annotation[which(as.vector(Experiment@Results@Annotation$Feature)==feat.names[i]),, drop=FALSE]
		if(nrow(subAnn)==0) next
		adduct.txt <- paste(apply(subAnn,1,function(x) {
			#x <- subAnn[1,]
			paste(as.matrix(x[c('Adduct','mIsoMass')]), collapse=' ')	
		}), collapse='; ')
		annid.txt <- paste(subAnn[,'AnnID'], collapse='; ')
		alignid.txt <- as.numeric(as.vector(subAnn[1,'AlignID']))
		iso.txt <- as.character(subAnn[1,'Isotope'])
		if(is.na(iso.txt)) iso.txt <- ''
		toMSMS.txt <- as.character(subAnn[1,'toMSMS'])
		if(is.na(toMSMS.txt)) toMSMS.txt <- ''
		
		fT[i, c('Adduct', 'AnnID', 'AlignID','Isotope','toMSMS')] <- c(adduct.txt, annid.txt, alignid.txt,iso.txt,toMSMS.txt)	
	}
	noAnItems <- match(as.vector(alignIDTab[,1]), fT[,1])
	fT[noAnItems, 'AlignID'] <- alignIDTab[,2]
	
	fT
}

create.matrix <- function(dim1,dim2) {
  x <- matrix()
  length(x) <- dim1*dim2
  dim(x) <- c(dim1,dim2)
  x
}

# # retEIC <- function(xcmSet, filePaths, index=NULL)
# {
	# nfiles <- length(filePaths)
	# scantimes <- list()
	
	  # if(nfiles == 1){
	    # #Single sample experiment
	    # if (file.exists(filePaths[1])) { 
	
	     # xraw <- xcmsRaw(filePaths[1],profstep=0)
	     # maxscans <- length(xraw@scantime)
	     # scantimes[[1]] <- xraw@scantime
	     # pdata <- as.data.frame(xcmSet@peaks) 
	
	     # EIC <- create.matrix(nrow(pdata),maxscans)
	     # #### IMPORTAR TAMBE AQUESTA FUNCIO!!
	     # EIC[,] <- getEIC4Peaks(xraw,pdata,maxscans)
	
	    # } else {
	      # stop('Raw data file:',filePaths[1],' not found. Cannot import data from XCMS. \n')
	    # }
	  # } else {
	    # #Multiple sample experiment
	    # gval <- groupval(xcmSet)
	    # na.flag <- 0
	    # maxscans <- 0	
	    # if (file.exists(filePaths[1])) { 
	      # xraw <- xcmsRaw(filePaths[1],profstep=0)
	      # maxscans <- length(xraw@scantime)
	    # } else {
	      # stop('Raw data file:',filePaths[1],' not found. Cannot import data from XCMS. \n');
	    # }
	
	    # #generate EIC Matrix
	    # EIC <- create.matrix(nrow(gval),maxscans)
	
	    # #loop over all samples
	    # for (f in 1:nfiles){
	      # idx.peaks <- which(index == f)
	      # if(length(idx.peaks) == 0) next
	      # if (file.exists(filePaths[f])) {
	        # #read sample
	        # xraw <- xcmsRaw(filePaths[f], profstep=0);
	        # maxscans.tmp <- length(xraw@scantime);
	        # scantimes[[f]] <- xraw@scantime
	        # if(maxscans.tmp > maxscans){
	          # #increase columns of EIC matrix
	          # EIC <- cbind(EIC,create.matrix(nrow(gval),maxscans.tmp - maxscans));
	          # maxscans <- maxscans.tmp;
	        # }
	        # pdata <- as.data.frame(xcmSet@peaks[gval[idx.peaks,f],,drop=FALSE]) # data for peaks from file f
	        # if(length(which(is.na(pdata[,1]))) > 0) na.flag <- 1
	        # EIC[idx.peaks,] <- getEIC4Peaks(xraw,pdata,maxscans)	
	      # } else {
	        # stop('Raw data file:',filePaths[f],' not found. Cannot import data from XCMS. \n')
	      # }
	    # }
	  # }
	  # list(scantimes=scantimes,EIC=EIC) 
# }

# getEIC4Peaks <- function(xraw,peaks,maxscans=length(xraw@scantime)){
  # if (!is.double(xraw@env$mz) || !is.double(xraw@env$intensity) || !is.integer(xraw@scanindex)) stop('mz/int not double.')
  # npeaks <- dim(peaks)[1]; 
  # scans  <- length(xraw@scantime);
  # eics <- matrix(NA,npeaks,maxscans);
  # for (p in 1:npeaks) {
    # timerange       <- c(peaks[p,"rtmin"],peaks[p,"rtmax"]);
    # tidx <- which((xraw@scantime >= timerange[1]) & (xraw@scantime <= timerange[2]));
    # if(length(tidx)>0){
      # scanrange <- range(tidx);
    # }else{
      # scanrange <- 1:scans;
    # }
    # massrange <- c(peaks[p,"mzmin"],peaks[p,"mzmax"]);
    # eic <- .Call("getEIC",xraw@env$mz,xraw@env$intensity,xraw@scanindex,as.double(massrange),
      # as.integer(scanrange),as.integer(length(xraw@scantime)), PACKAGE ='xcms' )$intensity;
    # eic[eic==0] <- NA;
    # eics[p,scanrange[1]:scanrange[2]] <- eic; 
  # }
# eics
# }


# show.xcmsSet <-  function(object) {
    # info.txt1 <- paste("A \"metaboSet\" object imported from an \"xcmsSet\" object with", nrow(object@phenoData), "samples\n\n", "Time range: ", paste(round(range(object@peaks[,"rt"]), 1), collapse = "-")," seconds (", paste(round(range(object@peaks[,"rt"])/60, 1), collapse = "-"), " minutes)\n", "Mass range:", paste(round(range(object@peaks[,"mz"], na.rm = TRUE), 4), collapse = "-"), "m/z\n", "Peaks:", nrow(object@peaks), "(about", round(nrow(object@peaks)/nrow(object@phenoData)), "per sample)\n", "Peak Groups:", nrow(object@groups), "\n", "Sample classes:", paste(levels(sampclass(object)), collapse = ", "), "\n\n")

	# info.txt2 <- NULL
	# info.txt3 <- NULL
	# info.txt4 <- NULL
    # if(.hasSlot(object, "mslevel")){
        # MSn <- mslevel(object)
        # if(is.null(MSn))
            # MSn <- 1
        # info.txt2 <- paste0("Peak picking was performed on MS", MSn, ".\n")
    # }
    # if(.hasSlot(object, "scanrange")){
        # if(!is.null(scanrange(object))){
            # info.txt3 <- paste("Scan range limited to ", scanrange(object)[1], "-", scanrange(object)[2], "\n")
        # }
    # }

    # if (length(object@profinfo)) {
        # for (i in seq(along = object@profinfo)) {
            # if (i != 1) cat("                  ")
            # info.txt4 <- paste("Profile settings: ", names(object@profinfo)[i], " = ", object@profinfo[[i]], "\n", sep = "")
        # }
    # }
	# info.txt <- paste(info.txt1, info.txt2, info.txt3, info.txt4)
	# info.txt
# }







