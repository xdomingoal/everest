\name{showAn}
\alias{showAn}
\title{Show annotation results from eveRest}
\description{Shows the detailed annotation results from a given feature of interest}
\usage{
showAn(Experiment, feature)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by evAnnotate.}
  \item{feature}{Character, the feature name for which detailed annotation results are queried.}
}
\details{
See eveRest vignette for more details. To open the vignette, execute the following code in R:
vignette("eveRestManual", package="everest")

}
\examples{
	#showAn(ex, 'M112T102') 
}

