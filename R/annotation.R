getPPM <- function(mass.error, mass) (mass.error*10^6)/round(mass)
getMassError <- function(ppm.error, mass) (round(mass)*ppm.error)/10^6

isotopeFilter <- function(MZ.vector, isotope.dis=1.003355, ppm.error=10)
{
	mzVect <- MZ.vector[order(MZ.vector, decreasing=F)]
	MassErr <- abs(abs(mzVect[-length(mzVect)] - mzVect[-1]) - isotope.dis)
	loc.prob <- sapply(1:length(MassErr), function(x) get.Prob(MassErr[x],0,round(mzVect[x])*(ppm.error/10^6)))
	if(length(which(loc.prob>0.6))==0) return(mzVect)
	mzVect[-(which(loc.prob>0.6)+1)]
}

isotopeFilterInt <- function(MZ.vector, INT.vector, isotope.dis=1.003355, ppm.error, featVect=NULL)
{
	#MZ.vector <- MZ.vect
	#INT.vector <- INT.vect
	
	## MZ check
	mzVect <- MZ.vector[order(MZ.vector, decreasing=F)]
	INT.vector <- INT.vector[order(MZ.vector, decreasing=F)]
	if(!is.null(featVect)) featVect <- featVect[order(MZ.vector, decreasing=F)]
	MassErr <- abs(abs(mzVect[-length(mzVect)] - mzVect[-1]) - isotope.dis)
	loc.prob <- sapply(1:length(MassErr), function(x) get.Prob(MassErr[x],0,round(mzVect[x])*(ppm.error/10^6)))
	if(length(which(loc.prob>0.6))==0){
		retMat <- matrix(c(INT.vector,rep(0,length(INT.vector)),rep(0,length(INT.vector)),mzVect), ncol=4)
		if(is.null(featVect)) return(retMat)
		if(!is.null(featVect)) return(cbind(retMat,featVect))
	}
	
	## Intensity Check
	dVect <- MZ.vector*0
	dVect[(which(loc.prob>0.6)+1)] <- 1
	#dVect <- c(0,1,1,0,0,0,1,0)
	ddVect <- c(diff(dVect),0)
	#ddVect[ddVect<0] <- 0
	dVect <- dVect + ddVect
	dVect[dVect==0] <- -1
	codontype <- dVect
	sqnc <- cumprod(codontype)*codontype
	sqnc.groups <- c(1,diff(sqnc)!=0)
	#cumsum(sqnc.groups)
	
	INT.mat <- cbind(INT.vector, cumsum(sqnc.groups), rep(0, length(INT.vector)), mzVect)
	
	#u <- sapply(unique(INT.mat[,2]), function(x){
	for(x in unique(INT.mat[,2]))
	{	
		y <- which(INT.mat[,2]==x)
		x.int <- INT.mat[y,1]
		xVar <- c(0,x.int[-length(x.int)]/x.int[-1])
		xVar[x.int<15] <- 2
		xVar[xVar>2] <- 2
		xVar[xVar!=2] <- 1
		xVar[1] <- 1
		INT.mat[y,3] <- xVar - 1
		INT.mat[y[(which(xVar==2)[1] - 1)],3] <- length(which(xVar==2))*-1
	}
	if(!is.null(featVect)) INT.mat <- cbind(INT.mat, featVect)
	
	#if(length(which(INT.mat[,3]==1))==0) return(mzVect)
	#mzVect[-which(INT.mat[,3]==1)]
	
	#mzVect[-(which(loc.prob>0.8)+1)]
	INT.mat
}



get.Prob <- function(Val, Mean, Sigma)
{
	z.par <- (Val-Mean)
	#mul.fun <- 1/sqrt(3.141593*2*Sigma^2)
	gaussian.function <- exp(-(z.par^2)/(2*Sigma^2))	
	gaussian.function
}


getAnnMat <- function(MZ.vect, INT.vect, annRul, ISOCountVect, max.mz.error, ppm.error, featVect){
	annMat <- matrix(0, nrow=0, ncol=8)
	colnames(annMat) <- c("AnnID","Feature","Mass","Adduct","IsoCount","Score","mIsoMass","minMz")
	idX <- 1
	for(mz.i in 1:length(MZ.vect))
	{
		## Removing low-abundant features for central evaluation.
		if(INT.vect[mz.i]<15) next	
		
		for(add.j in 1:nrow(annRul))
		{
			localAnnRul <- as.numeric(as.vector(annRul[,"massdiff"]))
			Mmass <- MZ.vect[mz.i]*(1/as.numeric(as.vector(annRul[add.j,"nmol"]))) - as.numeric(as.vector(annRul[add.j,"massdiff"]))
			annSp <- localAnnRul + rep(Mmass, length(localAnnRul)) * as.numeric(as.vector(annRul[,"nmol"]))
						
			mzInds <- sapply(MZ.vect, function(x) which.min(abs(annSp-x)))
			MassErr <- abs(annSp[mzInds] - MZ.vect)
			
			MZ.vectComp <- MZ.vect
			MZ.vectComp[which(MZ.vectComp<MZ.vect[mz.i])] <- MZ.vect[mz.i]
			loc.prob <- rep(0, length(MassErr))
			loc.prob[mz.i] <- 0
			if(any(MassErr[-mz.i]<=max.mz.error)) loc.prob <- sapply(1:length(MassErr), function(x) get.Prob(MassErr[x],0,round(MZ.vectComp[x])*(ppm.error/10^6)))
			loc.prob[loc.prob<0.6] <- 0
			loc.prob[mz.i] <- 0
			loc.prob[MassErr>max.mz.error] <- 0
			TotalProb <- prod(loc.prob[-which(loc.prob==0)])
			if(length(which(loc.prob!=0))==0) TotalProb <- 0
			#if(TotalProb!=0) stop()
			if(TotalProb!=0) {
	
				loc.prob.nz <- loc.prob[-which(loc.prob==0)]
				Nz <- length(loc.prob.nz)
				TotalProb <- 0
				TotalProb <- round(100*sum(loc.prob.nz^(1/Nz))/Nz,0)
							
				probMat <- c(idX, featVect[mz.i], round(MZ.vect[mz.i],6), as.vector(annRul[add.j, 'name']), ISOCountVect[mz.i], TotalProb, round(Mmass,6), round(MZ.vect[mz.i],6))
				
				for(i.x in which(loc.prob!=0)){					
					probMat <- rbind(probMat, c(idX, featVect[i.x], round(MZ.vect[i.x],6), as.vector(annRul[mzInds[i.x], 'name']), ISOCountVect[i.x], TotalProb, round(Mmass,6), round(MZ.vect[mz.i],6)))
				
				}
				
				alrAnn <- unlist(apply(probMat[,c(2,4)], 1, function(x) which(x[1]==annMat[,2] & x[2]==annMat[,4])))
				if(length(alrAnn)!=0){
						if(any(MZ.vect[mz.i] < annMat[alrAnn,'minMz'])) 
						{
							annMat <- annMat[-which(annMat[,1] %in% annMat[alrAnn,1][1]),]
							annMat <- rbind(annMat, probMat)
							idX <- idX + 1
						}
				}else{
					annMat <- rbind(annMat, probMat)
					idX <- idX + 1
				}
			}
		}
			
	}
	annMat
}

#@topAdducts: from the list of prioritized adducts, number of top adducts to be considered. 

annotateSP <- function(SpectralChain, annRul, ppm.error, featVect, max.mz.error, topAdducts, db.mz)
{
	
	#SpectralChain=as.vector(AL.list[4,]$Spectra)
	#featVect <- as.vector(strsplit(as.vector(AL.list[4,]$Features), ",")[[1]])	
			
	MZ.vect <- as.numeric(as.vector(sapply(strsplit(SpectralChain, " ")[[1]], function(x) strsplit(x, ",")[[1]][1]))) 
	INT.vect <- as.numeric(as.vector(sapply(strsplit(SpectralChain, " ")[[1]], function(x) strsplit(x, ",")[[1]][2]))) 
	INT.vect <- round(1000*(INT.vect/max(INT.vect)))

	if(length(which(INT.vect==0))!=0){
		MZ.vect <- MZ.vect[-which(INT.vect==0)]
		featVect <- featVect[-which(INT.vect==0)]		
		INT.vect <- INT.vect[-which(INT.vect==0)]
	}
	#MZ.vect
	annProb.list <- list()
	annP.vect <- vector()
	annP.table <- matrix(0, nrow=0, ncol=nrow(annRul)+6)
	colnames(annP.table) <- c("Feature","Mass","IsoCount", "M", "Add", "Prob",  as.vector(annRul[,"name"]))
	annC.vect <- vector()
	
	ISOmat <- isotopeFilterInt(MZ.vect, INT.vect, ppm.error=ppm.error, featVect=featVect)
	ISOmat <- cbind(ISOmat, rep(0, nrow(ISOmat)))
	colnames(ISOmat) <- c('INT.vector', 'isoid','isocount', 'mass', 'featVect','isoname')
	goOutMz <- which(ISOmat[,'isocount']==1)
	
	if(length(goOutMz)!=0)
	{
		MZ.vect <- as.numeric(as.vector(ISOmat[-goOutMz,'mass']))
		INT.vect <- as.numeric(as.vector(ISOmat[-goOutMz,'INT.vector']))
		if(!is.null(featVect)) featVect <- ISOmat[-goOutMz,'featVect']
		isoClm <- as.numeric(as.vector(ISOmat[,'isocount']))
		isoIndMat <- rbind(which(isoClm<0) + 1, c((which(isoClm<0) - 1)[-1], nrow(ISOmat)))
		isVect <- unlist(apply(isoIndMat, 2, function(x) {cmSm <- cumsum(isoClm[x[1]:x[2]])
			if(x[1]==x[2]) return(1)
			unique(cmSm)}))
		isVect.char <- sapply(isVect, function(x) paste('M+',x, sep=''))	
		ISOmat[goOutMz,'isoname'] <- isVect.char
	}else{
		featVect <- ISOmat[,'featVect']
		MZ.vect <- as.numeric(as.vector(ISOmat[,'mass']))
		INT.vect <- as.numeric(as.vector(ISOmat[,'INT.vector']))
	}
	ISOCountVect <- as.vector(-1*as.numeric(as.vector(ISOmat[which(ISOmat[,'isocount']<=0),'isocount'])))
	if(is.null(featVect)) featVect <- rep(0,length(MZ.vect))
	
	## Annotation Matrix
	annMat <- getAnnMat(MZ.vect, INT.vect, annRul, ISOCountVect, max.mz.error, ppm.error, featVect)
		
	## Annotation of the most intense m/z in MZ.vect
	
	# maxFeat <- which(annMat[,'Feature'] %in% featVect[which.max(INT.vect)])
	# if(length(maxFeat)==0) {
		# i.x <- which.max(INT.vect)
		# mostFreqAdd <- which.max(as.numeric(as.vector(annRul[,"freq"])))
		# Mmass <- round(MZ.vect[i.x],6) - as.vector(annRul[mostFreqAdd, 'massdiff'])
		# annMat <- rbind(annMat, cbind(idX, featVect[i.x], round(MZ.vect[i.x],6), as.vector(annRul[mostFreqAdd, 'name']), ISOCountVect[i.x], 'NoAdduct', round(Mmass,6), round(MZ.vect[i.x],6)))
	# }
	
	#withScore <- which(!is.na(annP.vect))

	## Adding isotopes if no annotation is given	
	if(nrow(annMat)==0 & any(as.numeric(as.vector(ISOmat[,'isocount']))<0)){
		#if(length(which(as.numeric(as.vector(ISOmat[,'isocount']))<0))>1) stop()
	
		isoElements <- which(as.numeric(as.vector(ISOmat[,'isocount']))<0)
		annIsotope <- as.data.frame(matrix('', nrow=0, ncol=10), stringsAsFactors=FALSE)
		colnames(annIsotope) <- c("AnnID", "Feature", "Mass", "Adduct", "IsoCount", "Score", "mIsoMass", "Isotope", "toMSMS", "featINT")

		for(i in 1:length(isoElements)){
			wIsoEl <- which(ISOmat[,'isoid'] %in% ISOmat[isoElements[i],'isoid'])
			#annIsotope <- rbind(annIsotope, matrix('', nrow=length(wIsoEl), ncol=10), stringsAsFactors=FALSE)

			add.txt <- ISOmat[wIsoEl,'isoname']
			add.txt[add.txt==0] <- 'M+0'
			isoVC <- as.numeric(as.vector(ISOmat[wIsoEl,'isocount']))
			isoFlag <- rep('', length(wIsoEl))
			isoFlag[isoVC>0] <- 'yes'
			isoVC[isoVC>0] <- 0
			isoVC <- abs(isoVC)
			isoVC[isoVC==0] <- ''
			
			isoLcMat <- cbind(rep(i, length(wIsoEl)), ISOmat[wIsoEl,c('featVect','mass')], add.txt, isoVC, rep('', length(wIsoEl)), rep('', length(wIsoEl)), isoFlag, rep('', length(wIsoEl)), ISOmat[wIsoEl,c('INT.vector')])
			colnames(isoLcMat) <- c("AnnID", "Feature", "Mass", "Adduct", "IsoCount", "Score", "mIsoMass", "Isotope", "toMSMS", "featINT")

			annIsotope <- rbind(annIsotope, isoLcMat)				
		}
		wIsoInds <- which(annIsotope[,'IsoCount']!='')
		toMSMSind <- wIsoInds[which.max(as.numeric(as.vector(annIsotope[wIsoInds,'featINT'])))]
		annIsotope$toMSMS <- factor(annIsotope$toMSMS, levels=c(levels(annIsotope$toMSMS), 'yes')) 
		annIsotope[toMSMSind,'toMSMS'] <- 'yes'
		annIsotope <- annIsotope[,-ncol(annIsotope)]
		return(annIsotope)
	} 
	
	if(nrow(annMat)==0) return("NoAnnotation")	

	annEnd <- annMat[,-ncol(annMat), drop=FALSE]
	rownames(annEnd) <- NULL

	## Adding isotopes (those that do not belong to an AnID should be also included here)
	if(length(goOutMz)!=0)
	{	
		annEnd.iso <- matrix(0, nrow=0, ncol=ncol(annEnd))
		colnames(annEnd.iso) <- colnames(annEnd)
		
		featIsotopes <- unique(annEnd[,"Feature"])[which(unique(annEnd[,"Feature"]) %in% ISOmat[which(as.numeric(as.vector(ISOmat[,3]))<0),'featVect'])]
		for(i in 1:length(featIsotopes))
		{
			lFeat <-  featIsotopes[i]
			for(j in unique(as.numeric(as.vector(annEnd[which(annEnd[,"Feature"] %in% lFeat),'AnnID']))))
			{
				locannEnd <- annEnd[which(annEnd[,'AnnID'] %in% j)[1],]
				isoID <- ISOmat[which(ISOmat[,'featVect'] %in% lFeat), 'isoid']
				locisomat <- ISOmat[which(ISOmat[,'isoid'] %in% isoID & ISOmat[,'isocount']>0 ), , drop=FALSE]
				locAnnEndMat <- t(matrix(rep(locannEnd, nrow(locisomat)), ncol=nrow(locisomat)))
				colnames(locAnnEndMat) <- colnames(annEnd)
				annEnd.iso <- rbind(annEnd.iso, cbind(locAnnEndMat[,'AnnID'],locisomat[,'featVect'], round(as.numeric(as.vector(locisomat[,'mass'])),5), locisomat[,'isoname'], 0, locAnnEndMat[,'Score'], locAnnEndMat[,'mIsoMass']))
				
			}	
		}
	rownames(annEnd.iso) <- NULL
	annEnd <- rbind(annEnd, annEnd.iso)
	}

	## Natural frequency sorting:
	if(is.null(annRul$freq)){
		warning('No frequency was listed in the annotation rules provided by the user. By doing so no adducts are not weighted, and more uncommon adducts such as (M+Na)+[-H4O2] have the same importance as [M+H]+ when listed. Take a look at the default everest annotation rules (execute: annRules) and add a column to your custom annotation rules.')
		AcAddFreqs <- rep(0, nrow(annEnd))
	}else{
		AddFreqs <- as.numeric(as.vector(sapply(as.vector(annEnd[,"Adduct"]), function(x) annRul[which(annRul[,"name"]==x),"freq"])))
		AddFreqs[is.na(	AddFreqs)] <- 0
		AcAddFreqs <- rep(0, length(AddFreqs))
		for(i in as.numeric(as.vector(unique(annEnd[,"AnnID"])))){
			NlnX <- which(as.numeric(as.vector(annEnd[,"AnnID"]))==i)
			AcScore <- sum(AddFreqs[NlnX])
			AcAddFreqs[NlnX] <- rep(AcScore,length(NlnX))		
		}
	}
	
	annEndSorted <- as.data.frame(annEnd[order(AcAddFreqs, decreasing=T),, drop=FALSE], stringsAsFactors=FALSE)
	
	## 1st Filter. Remove non-important adducts.
	
	isoFlag <- sapply(annEndSorted[,'Feature'], function(x) which(ISOmat[,'featVect'] %in% x))
	isoFlagv <- ISOmat[isoFlag,'isocount']
	isoFlagv[isoFlagv<0] <- 0
	annEndSorted <- cbind(annEndSorted, isoFlagv)

	#topAdducts <- 4
	tAdduct <- annRul[order(-annRul[,"freq"]),"name"][1:topAdducts]
	nAnnIDs <- unique(annEndSorted[,'AnnID'])
	topCandidateProvided <- FALSE
	
	removeInds <- vector()
	for(i in 1:length(nAnnIDs)){
		subTab <- annEndSorted[which(annEndSorted[,'AnnID']==nAnnIDs[i] & isoFlagv==0),]
		if(nrow(subTab)<=2 & topCandidateProvided) {
			removeInds <- c(removeInds, which(annEndSorted[,'AnnID']==nAnnIDs[i]))
			next
		}
		isTopAdduct <- sapply(subTab[,'Adduct'], function(x) x %in% tAdduct)
		if(any(isTopAdduct)) topCandidateProvided <- TRUE 
		if(!any(isTopAdduct) & topCandidateProvided){
			removeInds <- c(removeInds, which(annEndSorted[,'AnnID']==nAnnIDs[i]))
			next
		}

	}
	
 	if(length(removeInds)!=0) annEndSorted <- annEndSorted[-removeInds,]
	
	## 2nd Filter. Check by database feasibility (only if multiple m/z are reported)
	
	nAnnIDs <- unique(annEndSorted[,'AnnID'])
	removeInds <- vector()
	if(length(nAnnIDs)!=1 & !is.null(db.mz)){ 
		for(i in 2:length(nAnnIDs)){
			locIsoMass <- as.numeric(as.vector(annEndSorted$mIsoMass[which(annEndSorted[,'AnnID']==nAnnIDs[i])]))[1]
			loc.db.mz <- db.mz[which(db.mz>(locIsoMass-max.mz.error) & db.mz<(locIsoMass+max.mz.error))]
			loc.ppm <- sapply(loc.db.mz, function(x) getPPM(abs(locIsoMass-x),locIsoMass))
			inDB <- any(loc.ppm<ppm.error)
			if(inDB==FALSE) removeInds <- c(removeInds, which(annEndSorted[,'AnnID']==nAnnIDs[i]))
		}
	}
	
 	if(length(removeInds)!=0) annEndSorted <- annEndSorted[-removeInds,]	
	
	## Re-ordering data and flagging most plausible adduct
	nAnnIDs <- unique(annEndSorted[,'AnnID'])
	newAnID <- 1:length(nAnnIDs)
	newAnIDv <- annEndSorted[,'AnnID']
	for(i in 1:length(nAnnIDs)) newAnIDv <- gsub(nAnnIDs[i], newAnID[i], newAnIDv)
	annEndSorted[,'AnnID'] <- newAnIDv
	annEndSorted[,'isoFlagv'] <- gsub('1', 'yes', annEndSorted[,'isoFlagv'])
	annEndSorted[,'isoFlagv'] <- gsub('0', '', annEndSorted[,'isoFlagv'])
	annEndSorted <- cbind(annEndSorted, rep('',nrow(annEndSorted)), stringsAsFactors = FALSE)
	colnames(annEndSorted) <- c("AnnID", "Feature", "Mass", "Adduct", "IsoCount", "Score", "mIsoMass", "Isotope", "toMSMS")
	
	topAdductInd <- which(annEndSorted$AnnID==1 & annEndSorted$Isotope=='')
	topAdductValues <- annEndSorted[topAdductInd,'Adduct']
	AddFreqs <- as.numeric(as.vector(sapply(as.vector(topAdductValues), function(x) annRul[which(annRul[,"name"]==x),"freq"])))
	AddFreqs[is.na(	AddFreqs)] <- 0
	
	annEndSorted[topAdductInd[which.max(AddFreqs)],'toMSMS'] <- 'yes'
	
	annEndSorted
}

# formatUnAnnotated <- function(SpectralChain, featVect)
# {
	# MZ.vect <- as.numeric(as.vector(sapply(strsplit(SpectralChain, " ")[[1]], function(x) strsplit(x, ",")[[1]][1]))) 
	# INT.vect <- as.numeric(as.vector(sapply(strsplit(SpectralChain, " ")[[1]], function(x) strsplit(x, ",")[[1]][2]))) 
	# if(length(which(INT.vect==0))!=0) MZ.vect <- MZ.vect[-which(INT.vect==0)]
	# ISOmat <- isotopeFilterInt(MZ.vect, INT.vect, featVect=featVect)
	# goOutMz <- which(ISOmat[,3]==1)
	# if(length(goOutMz)!=0)
	# {
		# MZ.vect <- MZ.vect[-goOutMz]
		# if(!is.null(featVect)) featVect <- ISOmat[-goOutMz,4]
	# }else{
		# featVect <- ISOmat[,4]
	# }
	# ISOCountVect <- as.vector(-1*as.numeric(as.vector(ISOmat[which(ISOmat[,3]<=0),3])))
	# if(is.null(featVect)) featVect <- rep(0,length(MZ.vect))
	
	# OutMat <- cbind(featVect, round(MZ.vect,4), ISOCountVect, matrix("NoAnn", nrow=length(MZ.vect), ncol=5))
	# colnames(OutMat) <- c("Feature","Mass","IsoCount","M","Adduct","Prob","Mmean","SubGroup")

	# OutMat
# }

annotFeats <- function(AL.list, annRul, ppm.error, max.mz.error, db.object, fluidType, topAdductsN)
{
	#ppm.error <- 15
	#AL.list <- al.list
	#ppm.error 
	
	## DB/Fluid filtering
	if(!is.null(db.object)){
		if(is.null(fluidType)) db.mz <- as.numeric(as.vector(unlist(lapply(db.object@database, function(x) x$monisotopic_molecular_weight))))
		if(!is.null(fluidType)){
			fluidType <- tolower(fluidType)
			inFluid <- as.vector(unlist(lapply(db.object@database, function(x) (fluidType %in% tolower(x$biofluid_locations)))))
			db.mz <- as.numeric(as.vector(unlist(lapply(db.object@database[inFluid], function(x) x$monisotopic_molecular_weight)))) 
			infomsm <- paste('A total of', length(db.mz), 'metabolites found in fluid type:', fluidType, '\n')
			if(length(db.mz)==0) {
				warning('All metabolites in the database are going to be used since the selected biofluid had zero entries in the database. Make sure that the name of the selected biofluid matches exactly the one in the database.')
				db.mz <- as.numeric(as.vector(unlist(lapply(db.object@database, function(x) x$monisotopic_molecular_weight))))
			}else{
				cat(infomsm)
			}
		}			
	}else{
		db.mz <- NULL
	}
		
	pb <- txtProgressBar(min=1,max=nrow(AL.list), width=50, style=3)	
	annOList <- matrix(0, nrow=0, ncol=8)
	for(i in 1:nrow(AL.list))
	{	
		#which(as.numeric(as.vector(AL.list$AlignID))==344)
		
		featVect <- as.vector(strsplit(as.vector(AL.list[i,]$Features), ",")[[1]])
		locAnn <- annotateSP(SpectralChain=as.vector(AL.list[i,]$Spectra), annRul=annRul, ppm.error=ppm.error, featVect=featVect, max.mz.error=max.mz.error, topAdducts=4, db.mz=db.mz)
		
		if(unlist(locAnn)[1]=="NoAnnotation") next	
		
		if(is.null(dim(locAnn))) {
			locAnn <- t(as.data.frame(locAnn, stringsAsFactors=FALSE))
			rownames(locAnn) <- NULL
		}	
		
		alMat <- matrix(AL.list[i,]$AlignID, nrow=nrow(locAnn), ncol=1)
		colnames(alMat) <- "AlignID"
		annOList <- rbind(annOList,cbind(alMat,locAnn), stringsAsFactors=FALSE)		
		setTxtProgressBar(pb, i)
	}
	annOList <- as.data.frame(annOList)
	rownames(annOList) <- NULL
	annOList
}


# # db <- lapply(hmdb@database, function(x){
	
	# #x <- hmdb@database[[1]]
	# biofluid_locations <- unlist(x$biofluid_locations)
	# names(biofluid_locations) <- NULL
	# list(accession=x$accession, name=x$name, chemical_formula=x$chemical_formula, monisotopic_molecular_weight=x$monisotopic_molecular_weight, biofluid_locations=biofluid_locations, chemspider_id=x$chemspider_id, kegg_id=x$kegg_id, biocyc_id=x$biocyc_id, metlin_id=x$metlin_id, pubchem_compound_id=x$pubchem_compound_id, chebi_id=x$chebi_id)
	
# })

# hmdb@database <- db

# save(hmdb, file='hmdb.rda')

