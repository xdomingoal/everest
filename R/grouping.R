
groupIonsTable <- function(data.table, min.lin.rel, max.time.dist)
{
	featureTable <- data.table[,-1]	
	colnames(featureTable)[1:2] <- c("mz","RT")
	rownames(featureTable)=as.vector(data.table[,1])

	min.spectra.cor <- min.lin.rel	
	stopifnot(min.spectra.cor<1,min.spectra.cor>0)
	N.samples <- ncol(featureTable) - 2
	
	retention.time.vector <- as.vector(featureTable[,"RT"])
	ret.iterator <- as.matrix(1:length(retention.time.vector))
	
	## New Method:
	
	order.vector <- order(retention.time.vector)
	retention.time.vector.o <- retention.time.vector[order.vector]
	
	res.vector <- retention.time.vector.o - c(retention.time.vector.o[-1],retention.time.vector.o[length(retention.time.vector.o)])
	group.flags <- which(abs(res.vector)>max.time.dist)
	
	if(length(group.flags)==0)
    {
    	time.dist.clustlist <- list(order.vector[1:length(order.vector)])
    }else{
	    if (group.flags[1] != 1) group.flags <- c(0, group.flags)
	    if (group.flags[length(group.flags)] != length(order.vector))  group.flags <- c(group.flags, length(order.vector))
	    time.dist.clustlist <- sapply(1:(length(group.flags) - 1), function(i) order.vector[(group.flags[i] + 1):group.flags[(i + 1)]])
	}

	del.inds <- which(unlist(lapply(time.dist.clustlist, length))==1)
	if(length(del.inds)!=0) time.dist.clustlist <- time.dist.clustlist[-del.inds]

	###################
	#k <- 1
	pb <- txtProgressBar(min=0,max=length(time.dist.clustlist), width=50, style=3)
	global.aligned.factors <- list()
	for(k in 1:length(time.dist.clustlist))		
	{
		local.clust <- time.dist.clustlist[k][[1]]
		#class.vector <- global.class.vector[local.clust]
			
		time.dist <- as.matrix(dist(retention.time.vector[local.clust]), upper=T, diag=F)
				
		#local.spectra.matrix <- sapply(local.clust, function(x){
		#	convertMSPspectra(factors.list[[factors.assignment.matrix[x,1]]]$Spectra[[factors.assignment.matrix[x,2]]],max.mz)
		#})
				  
		cor.dist <- fastCor(t(featureTable[local.clust ,-c(1:2)]))
		#cor.dist <- suppressWarnings(cor(local.spectra.matrix))
		cor.dist[cor.dist<0] <- 0
		cor.dist <- 1 -  cor.dist

		cor.dist[abs(cor.dist)>(1-min.spectra.cor)] <- NA  ##Eliminar upper limits
		
		time.dist[abs(time.dist)>max.time.dist] <- NA
		Norm.factor <- max(time.dist, na.rm=TRUE)
		Norm.factor <- 1 # adding this
		time.dist <- time.dist/Norm.factor
		max.eu.dist <- sqrt((1-min.spectra.cor)^2+(max.time.dist/Norm.factor)^2)

		eu.dist <- sqrt(time.dist^2+cor.dist^2)
		eu.dist[eu.dist==0] <- NA
		eu.dist.vector <- which(eu.dist<max.eu.dist, arr.ind=T)
		#eu.dist <- cor.dist # adding this
	
		eu.dist.graph <- graph.data.frame(eu.dist.vector, directed = FALSE)
		eu.dist.clustlist <- split(unique(as.vector(eu.dist.vector)), clusters(eu.dist.graph)$membership)

		clustlist.unit.length <- which(as.vector(unlist(lapply(eu.dist.clustlist,length)))==1)
		if(length(clustlist.unit.length)!=0) eu.dist.clustlist <- eu.dist.clustlist[-clustlist.unit.length]

		aligned.factors <-  lapply(eu.dist.clustlist,function(clust){
			#clust <- eu.dist.clustlist[[4]]
			
			inter.distance.matrix <- eu.dist[clust,clust]
			
			#forbidden.combinations <- which((class.vector[clust] %*% t(class.vector[clust]))== class.vector[clust]^2, arr.ind=F)
			#inter.distance.matrix[forbidden.combinations] <- NA
			inter.distance.matrix[inter.distance.matrix==0] <- NA
			
			clusts <- comp.clusters(inter.distance.matrix, seq(1:nrow(inter.distance.matrix)))
			clusts <- lapply(clusts,function(x) {clust[as.vector(x$elements)]})
			clusts
		})
		
		aligned.factors <- unlist(aligned.factors, recursive = FALSE)
		aligned.factors <- lapply(aligned.factors, function(x) local.clust[x])

		global.aligned.factors <- c(global.aligned.factors,aligned.factors)

	setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
	}
	
	
	aligned.factors <- global.aligned.factors
				
	AlignID <- rep(0, nrow(featureTable))			
	for(i in 1:length(aligned.factors))
	{
		x <- aligned.factors[[i]] #[[1]]	
		AlignID[x] <- i
		#for(j in 1:length(loc[,1])) factors.list[[loc[j,1]]]$AlignID[loc[j,2]] <- i	
	}
	

	AL.list <- matrix(nrow=length(unique(AlignID))-1,ncol=4)
	for(AlID in unique(AlignID))
	{
		#AlID <- 19
		if(ncol(as.matrix(featureTable[which(AlignID==AlID),-c(1:2)]))==1) next
		sp.int <- as.numeric(as.vector(rowSums(as.matrix(featureTable[which(AlignID==AlID),-c(1:2)]))))
		sp.int <- round(normalize(sp.int)*1000)
		sp.mass.acc <- as.numeric(as.vector(featureTable[which(AlignID==AlID),'mz']))
		sp.text <- paste(apply(cbind(sp.mass.acc, sp.int), 1, function(x) paste(x, collapse=",")), collapse=" ")
		
		
		AL.list[AlID,1] <- AlID
		AL.list[AlID,2] <- round(mean(as.numeric(as.vector(featureTable[which(AlignID==AlID),2])))/60,3)
		AL.list[AlID,3] <- paste(rownames(featureTable[which(AlignID==AlID),]), collapse=",")
		AL.list[AlID,4] <- sp.text

		#if(length(sp.int)!=length(sp.mass.acc)) stop("Error: length of masses and intensities does not coincide")
		#sp.xcms <- rep(0,1000)
		#sp.xcms[sp.mass] <- sp.int/max(sp.int)
	
	}
	colnames(AL.list) <- c("AlignID", "RT", "Features", "Spectra")	
	AL.list <- as.data.frame(na.omit(AL.list))
	
	AL.list$Features
	
	
	AL.list

}

groupIonsXCMS <- function(xset, min.lin.rel, max.time.dist)
{
	#var.lin <- match.arg(var.lin, c("into", "maxo"), several.ok=FALSE)
	#xcmA <- groupval(xset, value=var.lin) #area
	
	## Add this:
	xcmA <- xcms::groupval(xset, value='maxo') #area
	xcmI <- xcms::groupval(xset, value='into') #area

	
	#RTMZ <- t(sapply(rownames(xcmA), function(x) as.numeric(as.vector(strsplit(as.character(x), "/")[[1]]))))
	#colnames(RTMZ) <- c("mz","RT")
	XCMt <- data.frame(xset@groups)
	rownames(XCMt)=xcms::groupnames(xset)
	featureTable <- cbind(XCMt$mzmed,XCMt$rtmed, xcmA)	
	featureTable.m <- cbind(XCMt$mzmed,XCMt$rtmed, xcmI)	
	colnames(featureTable)[1:2] <- c("mz","RT")
	
	min.spectra.cor <- min.lin.rel
	
	#cor.dist <- fastCor(t(featureTable[,-c(1:2)]))
	#cor.dist[cor.dist<0] <- 0
	
	
	stopifnot(min.spectra.cor<1,min.spectra.cor>0)
	
	N.samples <- ncol(featureTable) - 2
	
	
	retention.time.vector <- as.vector(featureTable[,"RT"])
	ret.iterator <- as.matrix(1:length(retention.time.vector))
	
	## New Method:
	
	order.vector <- order(retention.time.vector)
	retention.time.vector.o <- retention.time.vector[order.vector]
	
	res.vector <- retention.time.vector.o - c(retention.time.vector.o[-1],retention.time.vector.o[length(retention.time.vector.o)])
	group.flags <- which(abs(res.vector)>max.time.dist)
	
	#Fixed Bug #AL_290316:
	if(length(group.flags)==0)
    {
    	time.dist.clustlist <- list(order.vector[1:length(order.vector)])
    }else{
	    if (group.flags[1] != 1) group.flags <- c(0, group.flags)
	    if (group.flags[length(group.flags)] != length(order.vector))  group.flags <- c(group.flags, length(order.vector))
	    time.dist.clustlist <- sapply(1:(length(group.flags) - 1), function(i) order.vector[(group.flags[i] + 1):group.flags[(i + 1)]])
	}

	del.inds <- which(unlist(lapply(time.dist.clustlist, length))==1)
	if(length(del.inds)!=0) time.dist.clustlist <- time.dist.clustlist[-del.inds]

	###################
	#k <- 1
	pb <- txtProgressBar(min=0,max=length(time.dist.clustlist), width=50, style=3)
	global.aligned.factors <- list()
	for(k in 1:length(time.dist.clustlist))		
	{
		local.clust <- time.dist.clustlist[k][[1]]
		#class.vector <- global.class.vector[local.clust]
			
		time.dist <- as.matrix(dist(retention.time.vector[local.clust]), upper=T, diag=F)
				
		#local.spectra.matrix <- sapply(local.clust, function(x){
		#	convertMSPspectra(factors.list[[factors.assignment.matrix[x,1]]]$Spectra[[factors.assignment.matrix[x,2]]],max.mz)
		#})
				  
		cor.dist.1 <- fastCor(t(featureTable[local.clust ,-c(1:2)]))
		cor.dist.2 <- fastCor(t(featureTable.m[local.clust ,-c(1:2)]))
		cor.dist <- pmax(cor.dist.1,cor.dist.2)
		#cor.dist <- suppressWarnings(cor(local.spectra.matrix))
		cor.dist[cor.dist<0] <- 0
		cor.dist <- 1 -  cor.dist

		cor.dist[abs(cor.dist)>(1-min.spectra.cor)] <- NA  ##Eliminar upper limits
		
		time.dist[abs(time.dist)>max.time.dist] <- NA
		Norm.factor <- max(time.dist, na.rm=TRUE)
		Norm.factor <- 1 # adding this
		time.dist <- time.dist/Norm.factor
		max.eu.dist <- sqrt((1-min.spectra.cor)^2+(max.time.dist/Norm.factor)^2)

		eu.dist <- sqrt(time.dist^2+cor.dist^2)
		eu.dist[eu.dist==0] <- NA
		eu.dist.vector <- which(eu.dist<max.eu.dist, arr.ind=T)
		#eu.dist <- cor.dist # adding this
	
		eu.dist.graph <- graph.data.frame(eu.dist.vector, directed = FALSE)
		eu.dist.clustlist <- split(unique(as.vector(eu.dist.vector)), clusters(eu.dist.graph)$membership)

		clustlist.unit.length <- which(as.vector(unlist(lapply(eu.dist.clustlist,length)))==1)
		if(length(clustlist.unit.length)!=0) eu.dist.clustlist <- eu.dist.clustlist[-clustlist.unit.length]

		aligned.factors <-  lapply(eu.dist.clustlist,function(clust){
			#clust <- eu.dist.clustlist[[4]]
			
			inter.distance.matrix <- eu.dist[clust,clust]
			
			#forbidden.combinations <- which((class.vector[clust] %*% t(class.vector[clust]))== class.vector[clust]^2, arr.ind=F)
			#inter.distance.matrix[forbidden.combinations] <- NA
			inter.distance.matrix[inter.distance.matrix==0] <- NA
			
			clusts <- comp.clusters(inter.distance.matrix, seq(1:nrow(inter.distance.matrix)))
			clusts <- lapply(clusts,function(x) {clust[as.vector(x$elements)]})
			clusts
		})
		
		aligned.factors <- unlist(aligned.factors, recursive = FALSE)
		aligned.factors <- lapply(aligned.factors, function(x) local.clust[x])

		global.aligned.factors <- c(global.aligned.factors,aligned.factors)

	setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
	}
	
	
	aligned.factors <- global.aligned.factors
				
	AlignID <- rep(0, nrow(featureTable))			
	for(i in 1:length(aligned.factors))
	{
		x <- aligned.factors[[i]] #[[1]]	
		AlignID[x] <- i
		#for(j in 1:length(loc[,1])) factors.list[[loc[j,1]]]$AlignID[loc[j,2]] <- i	
	}
	

	AL.list <- matrix(nrow=length(unique(AlignID))-1,ncol=4)
	for(AlID in unique(AlignID))
	{
		#AlID <- 19
		if(ncol(as.matrix(featureTable[which(AlignID==AlID),-c(1:2)]))==1) next
		sp.int <- as.numeric(as.vector(rowSums(as.matrix(featureTable[which(AlignID==AlID),-c(1:2)]))))
		sp.int <- round(normalize(sp.int)*1000)
		sp.mass.acc <- as.numeric(as.vector(XCMt[which(AlignID==AlID),1]))
		sp.text <- paste(apply(cbind(sp.mass.acc, sp.int), 1, function(x) paste(x, collapse=",")), collapse=" ")
		
		
		AL.list[AlID,1] <- AlID
		AL.list[AlID,2] <- round(mean(as.numeric(as.vector(featureTable[which(AlignID==AlID),2])))/60,3)
		AL.list[AlID,3] <- paste(rownames(XCMt[which(AlignID==AlID),]), collapse=",")
		AL.list[AlID,4] <- sp.text

		#if(length(sp.int)!=length(sp.mass.acc)) stop("Error: length of masses and intensities does not coincide")
		#sp.xcms <- rep(0,1000)
		#sp.xcms[sp.mass] <- sp.int/max(sp.int)
	
	}
	colnames(AL.list) <- c("AlignID", "RT", "Features", "Spectra")	
	AL.list <- as.data.frame(na.omit(AL.list))
	
	AL.list$Features
	
	
	AL.list
}


comp.clusters <- function (hdist, classes) 
{
	
	#hdist <- inter.distance.matrix
	#classes <- seq(1:nrow(inter.distance.matrix))
	
    k <- nrow(hdist)
    n.class <- unique(classes)
    group.list <- list()
    it <- 1
    while (1) {
        min.list <- list()
        min.index <- 0
        min.previous <- NA
        for (j in 1:k) {
            #minn <- rep(NA, length(n.class))
           #w.minn <- minn
            # for (i in n.class) {
               	# #D <- hdist[j, which(classes == i)]
                # D <- hdist[j, i]
                # if (!all(is.na(D))) {
                  # minn[i] <- min(D, na.rm = T)
                  # w.minn[i] <- which(classes == i)[which.min(D)]
                # }
       		# }
       		minn <- as.vector(hdist[j,])
       		w.minn <- classes
       		w.minn[is.na(minn)] <- NA
       		
            #Su[[j]] <- minn
            #stop()
            group.distance <- sqrt(sum(minn^2, na.rm = T))/length(which(is.na(minn) == F))
            if (is.na(group.distance)) 
                group.distance <- 0
            min.list[[j]] <- list(dist = group.distance, elements = w.minn)
            if (which.min(c(group.distance, min.previous)) == 
                1 && group.distance != 0) {
                min.previous <- group.distance
                min.index <- j
            }
        }
        if (min.index == 0) 
            break
        g.elements <- as.vector(na.omit(min.list[[min.index]]$elements))
        g.elements <- c(g.elements, min.index)
        if (!any(is.null(unlist(group.list)))) {
            remove.i <- which((g.elements %in% unique(unlist(group.list))) == 
                T)
            if (length(remove.i) != 0) 
                g.elements <- g.elements[-remove.i]
        }
        group.list[[it]] <- list(elements = g.elements)
        it <- it + 1
        hdist[, g.elements] <- NA
        hdist[g.elements, ] <- NA
    #cat(it, ",")
    }
    group.list
}

