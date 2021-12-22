---
title: 'RHEIA: Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems'
tags:
  - Python
  - hydrogen-based systems
  - robust design optimization
  - uncertainty quantification
authors:
  - name: Diederik Coppitters
    orcid: 0000-0001-9480-2781
    affiliation: "1, 2"
  - name: Panagiotis Tsirikoglou
    orcid: xxx
    affiliation: 3
  - name: Ward De Paepe
    orcid: 0000-0001-5008-2946
    affiliation: 4
  - name: Konstantinos Kyprianidis
    orcid: 0000-0001-5008-2946
    affiliation: 5
  - name: Anestis Kalfas
    orcid: 0000-0001-5008-2946
    affiliation: 6
  - name: Francesco Contino
    orcid: 0000-0002-8341-4350
    affiliation: 1
affiliations:
 - name: Institute of Mechanics, Materials and Civil Engineering, Université catholique de Louvain
   index: 1
 - name: Fluid and Thermal Dynamics, Vrije Universiteit Brussel
   index: 2
 - name: Limmat Scientific AG
   index: 3
 - name: Thermal Engineering and Combustion Unit, University of Mons
   index: 4
 - name: Department of Automation in Energy and Environment, School of Business, Society and Engineering, Malardalen University
   index: 5
 - name: Department of Mechanical Engineering, Aristotle University of Thessaloniki
   index: 6
   

   
date: 21 December 2021
bibliography: paper.bib
---

# Summary

The massive deployment of intermittent renewable energy sources, such as solar and wind, 
requires energy storage over longer periods, to cover the energy demand at all times.
Among others, renewable energy storage in the form of hydrogen is gaining a lot of attention. 
The stored energy can be recovered in multiple sectors:
hydrogen can be converted back into electricity (power-to-power),
it can be used to produce low-carbon fuels (power-to-fuel),
and it can be used to fuel hydrogen vehicles (power-to-mobility).
The performance of these hydrogen-based energy systems is subject to uncertainties, 
e.g., the solar irradiance, the energy consumption of hydrogen-powered buses, the price of grid electricity.
Disregarding these uncertainties in the design process can result in a drastic mismatch between simulated and real-world performance, 
and thus lead to a *kill-by-randomness* of the system.
RHEIA provides a robust design optimization algorithm, which considers real-world uncertainties and yields design solutions that are least sensitive to these uncertainties in their performance, i.e., the robust design.
Moreover, models are included to evaluate the techno-economic and environmental performance of hydrogen in a power-to-fuel, power-to-power and power-to-mobility context.
When combined, RHEIA unlocks the robust designs for hydrogen-based energy systems.


# Statement of need

In design optimization studies of renewable energy systems, incorporating hydrogen is still an anomaly [@Eriksson2017]. 
Moreover, the optimization is often performed under the assumption of deterministic parameters (i.e., fixed, free from inherent variation).
Considering fixed values for model parameters in design optimization yields designs sensitive to real-world uncertainties
and results in a drastic mismatch between simulated and actual performances.
In fields different from energy systems, e.g., structural mechanics, aerospace and automobile applications, 
Robust Design Optimization (RDO) yielded robust designs by minimizing the variance on the performance [@orosz2020robust].
Consequently, alternative design solutions were proposed which provide a performance that is least sensitive to the random environment.
To ensure the computational tractability of RDO, surrogate modelling techniques achieve a promising computational efficiency
to quantify the mean and variance on the performance. Nevertheless, applications of such surrogate-assisted robust design optimization techniques are limited [@Chatterjee2017].
To fill these research gaps, RHEIA provides a multi-objective robust design optimization algorithm,
for which the uncertainty quantification is performed through a Polynomial Chaos Expansion (PCE) surrogate modelling technique.
In addition, RHEIA includes Python-based models for relevant valorization pathways of hydrogen: power-to-fuel, power-to-power and power-to-mobility.
The significant techno-economic and environmental uncertainties for these models are characterized based on scientific literature
and a method is included to gather climate data and demand data for the location of interest.
Finally, RHEIA allows connecting your own models to the RDO and uncertainty quantification algorithms as well.   

Simulation models exist which include the evaluation of hydrogen-based energy systems,
e.g., [INSEL](https://insel.eu/en/home_en.html), [EnergyPLAN](https://www.energyplan.eu/) and [TRNSYS](http://www.trnsys.com/).
Despite their extensive component model libraries, these models lack an optimization feature.
[HOMER Energy](https://www.homerenergy.com/products/pro/index.html) includes an optimization algorithm to design hybrid microgrids, including hydrogen system component models.
In Python, [Calliope](https://www.callio.pe/) [@pfenninger2018calliope] considers the optimization of multi-scale energy system models, where hydrogen is considered as a fuel in advanced gas turbines.
However, neither multi-objective problems, nor uncertainties during design optimization can be considered.

<!---
what about dakota? mention it, where are we different?
Separate modules for surrogate-assisted robust design optimization are present as well:
The DEAP package [@fortin2012deap] includes evolutionary optimization algorithms, while UQlab [@marelli2014uqlab] (Matlab) and ChaosPy [@feinberg2015chaospy] provide the PCE algorithm.  

A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
--->

The RDO framework has been applied on hydrogen-based 
energy systems, developed in Python: A directly-coupled 
photovoltaic-electrolyzer system [@Coppitters2019a] and a photovoltaic-battery-hydrogen system [@coppitters2020robust]. In addition, 
an Aspen Plus model of a power-to-ammonia model has been optimized with the framework[@verleysen2020can].
Other Aspen Plus models have been optimized as well through 
RHEIA: a micro gas turbine with a carbon capture plant [@giorgetti2019] and a 
micro gas turbine [@Paepe2019a]. Finally, uncertainty quantification has been
performed on an EnergyScope model [@rixhon2021role].


# Future work

The following improvements will be made in future versions of RHEIA:

- Including a sparse PCE algorithm, developed in our research group, to handle the curse-of-dimensionality for high-dimensional problems [@Abraham2017].
The sparse PCE algorithm has been proven effective in RDO for a photovoltaic-battery-hydrogen application [@coppitters2020robust]. 
To ensure a smooth inclusion of this sparse PCE algorithm in RHEIA, we built the ``pce`` module, 
as opposed to adopting an existing PCE package in Python, such as ChaosPy [@feinberg2015chaospy].
- Including optimization algorithm alternatives (e.g., Particle Swarm Optimization, Firefly Algorithm, Cuckoo Search), 
following the experience gained in our research group on using these algorithms in a surrogate-assisted RDO context [@Tsirikoglou2017].
- Adding additional models on hydrogen-based energy carrier production and utilization (e.g. ammonia, biomethane) in power-to-gas applications. 
- Including an adapted PCE to perform uncertainty quantification with imprecise probabilities, to distinguish between the importance of
epistemic and aleatory uncertainty on a parameter. We performed an RDO with imprecise probabilities on a photovoltaic-battery-heat pump system[@coppitters2021robust].

# Acknowledgements

The first author acknowledges the support of Fonds de la Recherche Scientifique - FNRS [35484777 FRIA-B2].

# References