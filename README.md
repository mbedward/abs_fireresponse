# R package: fireresponse

This is a prototype R package to assess and predict the status of fauna functional groups to 
alternative fire regimes. It is based on research by Victoria Reynolds and colleagues in the 
Applied Bushfire Science group, New South Wales Department of Planning and Environment (DPE).
The package is being developed jointly by DPE and the Centre for Environmental Risk Management
of Bushfires, University of Wollongong.

Fire regimes are defined in terms of three components: frequency (in the preceding fifty year 
period), severity of last fire, and time since last fire. The likely response of each functional
group over a range of values for each of these components is determined via an expert-elicitation
process in which experts estimate the lower, highest and most likely expected value of relative 
abundance given a component value (e.g. 5 years since last fire). To derive a group response to
a given fire regime, defined by discrete values for frequency, severity and time since fire, the 
bounded estimates at each component value are viewed as a triangular distribution from which a
vector of random samples is drawn. The product of the three component vectors is then calculated
to give a vector of random overall response values. Finally, a beta distribution is fitted to this
vector to represent the density of group overall response values to the given fire regime. Summary
statistics from the fitted distribution (e.g. median, mode or selected quantiles) can then be used
for spatial prediction. The distribution can also serve as the prior for Bayesian updating to
assess the likely effects of subsequent fires.
