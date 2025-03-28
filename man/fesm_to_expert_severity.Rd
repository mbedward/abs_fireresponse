% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster_prediction.R
\name{fesm_to_expert_severity}
\alias{fesm_to_expert_severity}
\title{Convert FESM integer fire severity values to scale used for expert-elicited
data}
\usage{
fesm_to_expert_severity(sev)
}
\arguments{
\item{sev}{Vector of integer values for FESM fire severity (0-5).}
}
\value{
A vector of integer severity values on the expert-elicited scale
  (0-8).
}
\description{
The FESM scale of fire severity comprises integers in the range 0-5, where 0
is unburnt; 2 is low; 3 is moderate; 4 is high; and 5 is extreme. This
function transforms these values to the 0-8 scale used in the expert-elicited
data for functional group responses, such that FESM values \code{(0,2,3,4,5)}
are mapped to expert values \code{(0,2,4,6,8)}.
}
\details{
Values of 1 (referred to as 'reserved class' in FESM documentation) do not
usually occur in FESM data and will be treated as unburnt (0) by this
function.
}
\examples{
# Convert the commonly used FESM values
fesm_to_expert_severity(c(0, 2, 3, 4, 5))

# Can also convert FESM value of 1 (not normally used)
fesm_to_expert_severity(1)

}
