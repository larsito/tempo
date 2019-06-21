\name{calc_sync}
\alias{calc_sync}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate synchrony indices
}
\description{
Function to calculate all or a subset of synchrony indices, including their
three term local quadrat variance (t3) version, and the decomposition of the
variance ratio index into Stotal, Strend, and Sdetrended.
}
\usage{
calc_sync(x, decompose = TRUE, indices = c("eta", "eta_t3", "eta_w",
  "eta_t3_w", "phi", "phi_t3", "varrat", "varrat_t3", "log_varrat",
  "log_varrat_t3"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{temporal community data in a data frame. Species as columns, years
as rows.}

\item{decompose}{If TRUE, the outputted synchrony indices will also contain
the total variance ratio (Stotal) and its decomposition into parts
attributable to long term trend (Strend) and actual synchrony (Sdetrended).
See function decompostr.}

\item{indices}{A character vector that specifies wich of the synchrony
indices should be calculated. By default, all available indices are
calculated. If only a subset is desired, the names of the functions to
calculate must be provided in the character vector (e.g. c("eta", "phi") to
calculate only standard versions of eta and phi).}
}
\value{
A data frame containing all indices that were specified to be
calculated by the indices and decompose arguments.
}
\details{
At the moment, the following indices of synchrony are available (in the order they are displayed in the function output, give that all indices haven been calculated):

\code{eta}: The mean correlation of abundances of each species with all other species.

\code{eta_w}: As eta, but the mean is weighted by the total abundances of species.

\code{eat_t3}, \code{eta_t3_w}: The t3 versions of eta and eta_w, respectively.

\code{phi}: Phi after Loreau and de Mazancourt (2008).

\code{phi_t3}: t3 modified version of Phi.

\code{varrat}: variance ratio after Schluter (1984), see also Leps ().

\code{varrat_t3}: t3 version of the variance ratio.

\code{log_varrat}: Log transformed variance ratio, see Leps et al. 2018.

\code{log_varrat_t3}: t3 version of the log variance ratio.

\code{syn_total}: The total synchrony, without decomposition (Stotal).

\code{syn_trend}: The component of total synchrony that is attributable to trend (Strend).

\code{syn_detrend}: The component of total synchrony that is attributable to year-to-year fluctuations (Sdetrend).
}
\references{
Blüthgen, N., N. K. Simons, K. Jung, D. Prati, S. C. Renner, S.
  Boch, M. Fischer, N. Holzel, V. H. Klaus, T. Kleinebecker, M. Tschapka, W.
  W. Weisser, and M. M. Gossner. 2016. Land use imperils plant and animal
  community stability through changes in asynchrony rather than diversity.
  Nat Commun 7:10697.

Gross, K., B. J. Cardinale, J. W. Fox, A. Gonzalez, M. Loreau, H.
  W. Polley, P. B. Reich, and J. van Ruijven. 2014. Species richness and the
  temporal stability of biomass production: a new analysis of recent
  biodiversity experiments. Am Nat 183:1-12.

Leps, J., M. Majekova, A. Vitova, J. Dolezal, and F. de Bello.
  2018. Stabilizing effects in temporal fluctuations: management, traits, and
  species richness in high-diversity communities. ECOLOGY 99:360-371.
  
Loreau, M., and C. de Mazancourt. 2008. Species Synchrony and Its
  Drivers: Neutral and Nonneutral Community Dynamics in Fluctuating
  Environments. The American Naturalist 172:E48-E66.

Schluter, D. 1984. A Variance Test for Detecting Species
  Associations, with Some Example Applications. ECOLOGY 65:998-998.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
