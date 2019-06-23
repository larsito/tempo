\name{syngenr}
\alias{syngenr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate fluctuations in species abundances across time
}
\description{
\code{syngenr} simulates the fluctuations of abundances of species in a community through time. Arguments can be set so that different fluctuation patterns can be simulated, i.e. synchrony, antisynchrony, and random fluctations. In addition, directional trends can be build into the species abundance fluctuations.

\code{syngenr} comes with predefined default values for all its arguments, so that running it without any arguments specified (i.e. \code{syngenr()}) will also return a result. 
}
\usage{
syngenr(years = 100, n_sp = 16, max_rel_abu = 0.6, tot_abu = 300, power = 1.8, 
switch_env = c("on", "off"), mean_env_resp = 1, sd_env_resp = 1, bimodal_env = FALSE, 
comp = FALSE, switch_trend = c("on", "off"), mean_trend_resp = 1, sd_trend_resp = 1, 
bimodal_trend = FALSE, bound_pos = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{years}{the length of the timeseries in years.}

\item{n_sp}{number of species in the community.}

\item{tot_abu}{total abundance of the community, representing e.g. the number of
individuals or the total amount of biomass.}

\item{power}{the slope of the relationship between log(mean) and log(variance) of
the abundances of the species.}

\item{switch_env}{character vector of length 1, either "on" or "off". Defines if the
species abundances repond to a hypothetical environmental cue.}

\item{mean_env_resp}{the mean of the normal distribution from which each of the
species responses to the environemtnal cue is drawn.}

\item{sd_env_resp}{the standard deviation around the mean of the normal distribution
from which each of the species responses to the environemtnal cue is drawn.}

\item{bimodal_env}{Logical. Making the reponse to the environmental cue either
uniform among the species (if FALSE), or making the majority of half of the
species respond positively, and the other half negatively.}

\item{comp}{Logical. If TRUE, species exhibit compensatory dynamics, i.e. the gain
of abundance in a species from one year to the next, is compensated by the
loss of abundance in another species, where the latter has a similar mean
abundance value.}

\item{switch_trend}{Either "on" or "off", defining if there is a general monotonic
trend in abundances of species across the timeseries.}

\item{mean_trend_resp}{the mean of the normal distribution from which each of the
species responses to the longterm trend is drawn. Default is 1.}

\item{sd_trend_resp}{the standard deviation around the mean of the normal
distribution from which each of the species responses to the longterm trend is
drawn. Default is 1.}

\item{bimodal_trend}{Logical. If TRUE, the majority of half of the species exhibit a
positive long term trend of abundances, and the other half expresses a
negative long term trend. If FALSE, most of the species exhibit a uniform long
term trend in abundances, depending on the value set for mean_trend and
mean_sd.}

\item{bound_pos}{Logical. If true, abundance values that are simulated to be
negative, will be set to zero.}
}
\details{

}
\value{
The output is a list with four elements

\code{$time_species_matrix}: the simulated temporal community data where species are columns and years are rows.

\code{$param_years}: values for the evnironmental cue and the trend throughout the years. Note that these contain values even if the environment or trend are
switched off.

\code{$param_species}: a data frame containing the responses to the environment and the long term trend, as well as the mean abundance and its standard deviation for each species in the community.

\code{$param_general}: a data frame with only one row, containing all the parameter settings from the function call.
}
\references{
Leps, J., Götzenberger, L., Valencia E., de Bello, Francesco. 2019. Accounting for long-term   
  directional trends on year-to-year synchrony in species fluctuations. Ecography. 
}
\author{
Lars Götzenberger
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Using syngenr with it's default values will simulate community of 16 species and 
# 100 years, with a signature of synchrony and a common positive longterm trend of species. 
simcom_sync_trend <- syngenr()

# The actual community is saved in the first list element 
simcom_sync_trend[[1]]

# Community with 8 species, 40 years, anti-synchrony, and not long-term trend
simcom_async_notrend <- syngenr(years = 40, n_sp = 8, bimodal_env = TRUE, comp = TRUE,
switch_trend = "off")
simcom_async_notrend[[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
