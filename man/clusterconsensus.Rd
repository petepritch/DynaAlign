% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterbreak.R
\name{clusterconsensus}
\alias{clusterconsensus}
\title{Generate consensus sequence}
\usage{
clusterconsensus(df)
}
\arguments{
\item{df}{Input clusterbreak output df}
}
\value{
df with first column being unique cluster id and second column being the corresponding consensus sequence
}
\description{
Generate consensus sequence
}
\examples{
library(DynaAlign)
test <- h3n2sample \%>\% dplyr::group_by(clade) \%>\% dplyr::sample_frac(.4)
test <- test \%>\% dplyr::distinct(sequence,.keep_all=T)
out.df <- clusterbreak(h3n2sample$sequence,size_max = 800,thresh_p=.8,sim_fn=function(x) similarityMH(x,k=4,n_hash=500))
clusterconsensus(out.df$clustered_seq)
}
