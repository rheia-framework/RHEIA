.. _lab:functionalities:

A short description
===================

RHEIA is a computational framework that facilitates multi-objective deterministic and robust design optimization of engineering models, 
as well as uncertainty quantification.
The deterministic design optimization process (RHEIA uses the well established Nondominated Sorting Genetic Algorithm, NSGA-II) considers 
deterministic values for any model parameter (i.e. parameters which are not directly controlled by the designer) in order to define an optimal model configuration. 
In robust design optimization, the model parameters can be considered uncertain (e.g. the uncertainty on the future evolution of the grid electricity price). 
As these model parameters are uncertain, the quantities of interest (i.e. the model outputs of interest, e.g. the levelized cost of electricity) are uncertain as well. 
Consequently, the objectives of the robust design optimization are the optimization of the mean and the minimization of the standard deviation of one or more outputs of interest. 
While optimizing the mean corresponds to optimizing the stochastic performance, minimizing the standard deviation corresponds to minimizing the sensitivity of the considered outputs of interest in respect to the uncertain model parameters. 
The uttermost goal of this process is to derive a full spectrum of model configurations ranging from the optimal stochastic performance to the least sensitive (i.e. most robust) one. 
Comparisons of the optimal model configurations as they are obtained from both deterministic and robust design optimization processes can highlight the effect of model parameters uncertainties to the desired model performance. 
From the algorithmic standpoint, to determine the mean and standard deviation on the  model output of interest, each set of uncertain model parameters (i.e. design sample) is subject to an uncertainty propagation process. 
Details on the robust design optimization procedure are provided in :ref:`lab:rdoprocedure`.
In RHEIA the uncertainty quantification is facilitated through the Polynomial Chaos Expansion (PCE).
In addition to the mean and standard deviation of any model output considered, PCE can provide a full characterization of the outputs' probability density function,
cumulative distribution function. Ranking of the model parameters, in terms of their contribution to the variance of the model outputs is also feasible, using Sobol's indices metric.   

As regards the end-use of RHEIA framework, let us first clarify that RHEIA can be coupled to any well-defined, computational model (a Python-wrapper is provided to connnect any open- or closed-source model). However, a clear focus towards hydrogen-based energy system models is adopted. This is due to our belief that towards the transformation of the energy sector to a decarbonised setting, deterministic and robust design optimization has a key role to provide optimal configurations and enable robust transitions. 
At the current situation and despite its potential as an energy carrier, hydrogen is scarsely integrated in design optimization studies of hybrid renewable energy systems.
Moreover, when hydrogen is considered, mainly deterministic model parameters (i.e. perfectly known and free from inherent variation) are assumed, despite the uncertainty
that affects the performance of such systems (e.g. market values, climate conditions). Therefore, RHEIA combines robust design optimization with hydrogen-based energy systems, to unlock robust designs which are least-sensitive to the stochastic environment. 

The hydrogen-based energy system models featured in RHEIA correspond to the main hydrogen valorization pathways: power-to-power, power-to-mobility and power-to-fuel.
For each model, a specific set of parameters are defined with uncertainty, which are considered during robust design optimization and uncertainty quantification. 
The techno-economic and environmental uncertainties are adopted from scientific literature.
The models depend on the climate data and/or energy demand data. To anticipate this necessity, RHEIA provides a method to gather climate data and energy demand data for the location of interest (e.g. determining the optimized mean and standard deviation of the levelized cost of hydrogen for a photovoltaic-electrolyzer system in Brussels, Belgium). 
Evidently, the climate and demand database can be extended by your own climate and energy demand data. The full set of features provided aims to identify RHEIA as a decision support tool for optimal, yet robust, pathways. 