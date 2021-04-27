.. _lab:functionalities:

A short description
===================

RHEIA provides a tool to perform multi-objective deterministic design optimization, robust design optimization and uncertainty quantification.
The deterministic design optimization algorithm (Nondominated Sorting Genetic Algorithm, NSGA-II) considers deterministic values for the model parameters
(i.e. parameters which are not controlled by the designer). In robust design optimization, the model parameters can be considered uncertain
(e.g. the uncertainty on the future evolution of the grid electricity price). As these model parameters are uncertain, the quantities of interest
( i.e. the model outputs of interest, e.g. the levelized cost of electricity) are uncertain as well. Consequently, the objectives of the robust design optimization are the 
optimization of the mean and the minimization of the standard deviation of one or more quantities of interest. While optimizing the mean corresponds to
optimizing the stochastic performance, minimizing the standard deviation corresponds to minimizing the sensitivity of the quantity of interest to the random environment.
Hence, in this approach, the robust design can be characterized (i.e. the design which achieves the lowest standard deviation on the quantity of interest).
The advantages of this robust design can then be compared with the optimized deterministic designs and designs with an optimized mean on the quantity of interest. 
To determine the mean and standard deviation on the quantity of interest for each design sample evaluated during the robust design optimization procedure, an uncertainty quantification is performed (Polynomial Chaos Expansion).
In addition to the mean and standard deviation, the Polynomial Chaos Expansion provides the probability density function,
cumulative distribution function and the Sobol' indices on the quantity of interest.  

While the optimization and uncertainty quantification can be performed on any system model, RHEIA includes hydrogen-based energy system models.
Despite its potential as an energy carrier, hydrogen is scarsely integrated in design optimization studies of hybrid renewable energy systems.
Moreover, when hydrogen is considered, mainly deterministic model parameters (i.e. perfectly known and free from inherent variation) are assumed, despite the uncertainty
that affects the performance of such systems (e.g. market values, climate conditions).
Therefore, RHEIA combines robust design optimization with hydrogen-based energy systems, to unlock robust design which are least-sensitive to the random environment.
Moreover, the uncertainty quantification enables to define the Sobol' indices for an optimized design,
which illustrate the uncertainties that dominate the variance of the quantity of interest.
The hydrogen-based energy system models encollapse the main hydrogen valorization pathways: power-to-power, power-to-mobility and power-to-fuel.
For each model, a specific set of parameters are defined with uncertainty, which are considered during robust design optimization and uncertainty quantification. 
The techno-economic and environmental uncertainties are adopted from scientific literature.
The models depend on the climate data and/or energy demand data. The framework provides a method to gather climate data and energy demand data for the location of interest. This enables the user to use the framework as a decision support tool, 
e.g. determining the optimized mean and standard deviation of the levelized cost of hydrogen for a photovoltaic-electrolyzer system in Brussels, Belgium. 
Evidently, the climate and demand database can be extended by your own climate and energy demand data.

In addition to the hydrogen-based system models, the framework allows to connect your own model. 
Through the Python wrapper, RHEIA enables to connect Python-based, open-source and closed-source system models.
