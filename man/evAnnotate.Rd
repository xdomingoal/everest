\name{evAnnotate}
\alias{evAnnotate}
\title{Annotation by Everest}
\description{Annotation of features detected by XCMS or a list of features detected by other software}
\usage{
evAnnotate(xcmsSet=NULL, data.table=NULL, 
ion.mode=c('pos','neg'), min.correlation=0.7, 
max.time.dist=1, ppm.error=20, max.mz.error=0.05,
annRules=everestRules, db.object=hmdb, fluidType=NULL,
topAdductsN=4, import.profiles=FALSE,
alternative.data.path=NULL)
}
\arguments{
  \item{xcmsSet}{A 'xcmsSet' S4 object containing the experiment data previously created by XCMS.}
  \item{data.table}{A data table containing the list of features (see Details)}
  \item{ion.mode}{Character. Wether positive 'pos' or negative 'neg' type of study.} 
  \item{min.correlation}{The mininmum abundance correlation among samples for which pairs of features are considered to be from the same compound.} 
  \item{max.time.dist}{The maximum retention time distance (in seconds) for which pairs of features might be considered to be from the same compound.} 
  \item{ppm.error}{The m/z error un parts per milion (ppm) in which annotations are performed} 
  \item{max.mz.error}{Despite the ppm.error, the maximum absolute error (in Da) for adduct annotation. } 
  \item{annRules}{Annotation Rules (see Details)} 
  \item{db.object}{Database used for result-filtering. Use NULL if no filtering is to be applied.} 
  \item{fluidType}{If a database is used for filtering (db.object), select type of fluid.} 
  \item{topAdductsN}{Integer, cut-off number of prioritized adduct that remain as top-candidate adducts.} 
  \item{import.profiles}{If TRUE, EIC are imported from the XCMS object (currently not implemented).} 
  \item{alternative.data.path}{Alternative path when the data path is not the same as in the XCMS object (currently not implemented).} 
}
\details{
See eveRest vignette for more details. To open the vignette, execute the following code in R:
vignette("eveRestManual", package="everest")
}
\value{
The function returns an updated S4 'MetaboSet' class, where the LC-MS peaks have been annotated.
}
\references{
[1] Xavier Domingo-Almenara, et al., eveRest: A Computational Tool for Metabolite Annotation in LC/MS-Based Metabolomics. Anal. Chem (2018) Submitted.
}
\author{
Xavier Domingo-Almenara. xdomingo@scripps.edu
}

\seealso{
\code{\link{showAn}}, \code{\link{annoTable}}
}
\examples{

## 1a. Annotate compounds (from XCMS):

	#ex <- evAnnotate(xcmsSet=xset3, data.table=NULL, 
	#ion.mode='pos', min.correlation=0.6, max.time.dist=1, 
	#ppm.error=20, max.mz.error=0.05, annRules=annRules, 
	#import.profiles=FALSE, alternative.data.path=NULL)

## 1b. Annotate compounds (from Data Table):

	#ex <- evAnnotate(data.table=sample.tab, ion.mode='pos', 
	#min.correlation=0.6, max.time.dist=1, ppm.error=20, 
	#max.mz.error=0.05, annRules=annRules)

	#showAn(ex, 'M112T102') 
	#anTab <- annoTable(ex) 

## 2. Identify compounds:

	#ex <- evIdentify(ex, db=hmdb, ppm.error=5)

	#showId(ex, 'M112T102') 	
	#idTab <- idenTable(ex) 

}


