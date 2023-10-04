.. figure:: images/logo.svg
   :width: 80%
   :align: center

.. toctree::
   :maxdepth: 3
   :numbered:
   :hidden:
   
   functionalities
   connecting your own model
   installation
   tutorial
   auto_examples/index
   stochastic design space
   optimization
   uncertainty quantification
   energy system models
   methods
   API
   contribution
   bibliography
   glossary
   whatsnew

latest stable version: **v1.1.11**

Introduction
============

The Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems (RHEIA) framework provides 
multi-objective optimization (deterministic and stochastic) and uncertainty quantification algorithms. 
These algorithms can be applied on hydrogen-based energy systems, which are included in RHEIA.
In addition, RHEIA allows to connect your own models to the algorithms as well.

A brief overview on the features of RHEIA is provided in :ref:`lab:functionalities`,
followed by a detailed illustration on how to connect your own model (Python based, open source or closed source) in :ref:`lab:connectingyourownmodel`.
If these features comply with your need, the installation procedure and package dependencies are illustrated in :ref:`installationlabel`. 
As a first step, the :ref:`lab:tutorial` provides an initiation of using the framework on a hydrogen-based energy system. 
Additional examples are illustrated in :ref:`lab:examples`. 
The models are characterized by a design space, which defines the design variables, and a stochastic space when parameter uncertainty is considered.
Those spaces are defined in two files, which are elaborated in :ref:`lab:stochasticdesignspace`.
The guides to perform the deterministic design optimization, robust design optimization and uncertainty quantification are present in :ref:`lab:optimization` and 
:ref:`lab:uncertaintyquantification`, respectively.
The hydrogen-based energy system models are described in :ref:`lab:energysystemmodels`, including the uncertainty characterization of the techno-economic and environmental parameters, a script to evaluate the performance of specific designs and a method to get the climate data and demand data for your location of interest. The documentation concludes with brief details on the optimization and uncertainty quantification algorithms (:ref:`lab:methods`), the :ref:`lab:APIref` and the details on how to contribute to the framework (:ref:`lab:contribution`).

Citing RHEIA
============

Please use the following to cite the framework::

   Coppitters et al., (2022). RHEIA: Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems. Journal of Open Source Software, 7(75), 4370, https://doi.org/10.21105/joss.04370

Bibtex entry::

   @article{Coppitters2022,
     doi = {10.21105/joss.04370},
     url = {https://doi.org/10.21105/joss.04370},
     year = {2022},
     publisher = {The Open Journal},
     volume = {7},
     number = {75},
     pages = {4370},
     author = {Diederik Coppitters and Panagiotis Tsirikoglou and Ward De Paepe and Konstantinos Kyprianidis and Anestis Kalfas and Francesco Contino},
     title = {{RHEIA: Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems}},
     journal = {Journal of Open Source Software}
   }


Support
=======

Contact us for questions and feedback:

rheia.framework@gmail.com

License
=======

The project is licensed under `MIT license <https://github.com/rheia-framework/RHEIA/blob/main/LICENSE>`_. 


Indices and tables
==================

:ref:`genindex`
:ref:`modindex`
:ref:`search`
