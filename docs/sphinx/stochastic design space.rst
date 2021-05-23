.. _lab:stochasticdesignspace:

Define the design and stochastic spaces
=======================================

In this section, the characterization of the model parameters and design variables is described.
The design variables are variables that are controllable by the designer within a certain range. 
This range should be provided, to shape the search space for the optimizer. 
A model parameter corresponds to a parameter that usually cannot be controlled by the designer (e.g. the cost of a photovoltaic panel), 
or the decision on this parameter value is fixed (e.g. a fixed amount of photovoltaic panels on the roof). 
Such a parameter can be considered deterministic or uncertain.
The design variables and model parameters are characterized in two files: :file:`design_space` and :file:`stochastic_space`.
These files are present in the case folder for the specific case considered (e.g. :file:`CASES\\H2_FUEL` for the `H2_FUEL` case).
In :file:`design_space`, the deterministic values for the model parameters and the range for the design variables are provided.
In :file:`stochastic_space`, the uncertainty is allocated to the specific model parameters and design variables.
When a deterministic design optimization is performed, only the :file:`design_space` file is required. 
In the other cases, i.e. uncertainty quantification and robust design optimization, both files are required.

.. _lab:ssdesignspace:

The design_space file
---------------------

In the :file:`design_space` file, the design variables and the model parameters which need a quantification in your model are defined. 
When performing uncertainty quantification, the :file:`design_space` file consists only of model parameters.
In the case of deterministic design optimization or robust design optimization, the :file:`design_space` file requires design variables. 
Additionally, if some model parameters require a quantification outside the model, 
the :file:`design_space` file includes both :ref:`lab:ssdesignvariables` and :ref:`lab:ssmodelparameters`.

.. _lab:ssdesignvariables:

Characterizing the design variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
To define a design variable, the set-up in the :file:`design_space` file is as follows::

	name feature_type lb ub

where:

- name: name of the variable;
- feature_type: a design variable is indicated with `var`;
- lb: lower bound value for the design variable;
- ub: upper bound value for the design variable. 

An example for a configured design variable `des_var_1` with a range between 10 and 50 is::

    des_var_1 var 10 50


.. _lab:ssmodelparameters:

Characterizing the model parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this file :file:`design_space`, the deterministic value (or mean value when the parameter is considered uncertain) is provided.
The configuration of a model parameter is similar to the configuration of a design variable::

    name feature_type value

where:

- name: name of the variable;
- feature_type: a parameter is indicated with `par`;
- value: deterministic value (or mean value when the parameter is stochastic).

An example of a configured model parameter `par_1` with a mean value of 0.03 is::

	par_1 par 0.03

.. _lab:ssexampleds:

Example of design_space
^^^^^^^^^^^^^^^^^^^^^^^
An example of a configured :file:`design_space` file, which consists of 3 model parameters (par_1, par_2 and par_3) and 2 design variables (design_var_1 and design_var_2), is presented::

	design_var_1 var 1 3
	design_var_2 var 10 100
	par_1        par 4
	par_2        par 2.5
	par_3        par 175

.. _lab:ssstochastic_space:

The stochastic_space file
-------------------------

The uncertainty on the design variables and model parameters can be allocated through the file :file:`stochastic_space`. 
This file is required when performing robust design optimization and uncertainty quantification, where several parameters are subjected to uncertainty. 
For every design variable and model parameter defined in :file:`design_space`, an uncertainty can be defined.

Characterizing the uncertainties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Defining the uncertainty can be done as follows::

	name abs_rel distribution deviation

where:

	- name: name of the parameter or variable, equal to the name of the parameter or variable in :file:`design_space`;
	- abs_rel: absolute or relative uncertainty to the mean, defined with `absolute` or `relative`, respectively;
	- distribution: The distribution of the uncertainty;
	- deviation: uncertainty on the mean.

The following distributions are available:

- Uniform
- Gaussian

The meaning of deviation at the end of the line depends on the distribution. When a Uniform distribution is considered,
the deviation refers to the absolute (or relative) difference between the upper bound of the Uniform distribution and the mean: for :math:`\mathcal{U}(a,b)`, :math:`deviation = (b-a)/2`). 
When a Gaussian distribution is considered, the value corresponds to the standard deviation: :math:`\mathcal{N}(mean,deviation)`.
Keep always in mind that the mean value of the design variables is the deterministic value assigned by the optimizer in each iteration. In the case of model fixed parameter the mean value corresponds to the fixed value as it is assigned in the model definition.

An example of a configured uncertain parameter `par_2`, characterized by a Uniform distribution with a :math:`\pm 1` deviation from the mean value::

	par_2 absolute Uniform 1

Note that it is not required to allocate an uncertainty to every design variable and model parameter defined in :file:`design_space`.
In other words, when a parameter (or variable) is defined in :file:`design_space`, but not in :file:`stochastic_space`, the parameter (or variable) is considered deterministic. 
Moreover, it is not necessary to keep the same order of appearance of parameters and variables in :file:`design_space` :file:`stochastic_space` files.

Example of stochastic_space
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In summary, a :file:`stochastic_space` file corresponding to the illustrative :file:`design_space` example file in :ref:`lab:ssexampleds` looks like this::

	par_1        relative Gaussian 0.5
	par_2        absolute Uniform  1
	design_var_2 relative Uniform  0.1

Where the model parameter `par_3` and design variable `design_var_1` are considered deterministic, 
`par_1` is characterized by a Gaussian distribution with a 
relative standard deviation of 0.5 (i.e. :math:`\mathcal{N}(4,2)`),    
`par_2` is characterized by a Uniform distribution with an 
absolute deviation of 1 (i.e. :math:`\mathcal{U}(1.5,3.5)`) and    
`design_var_2` is characterized by a Uniform distribution with a 
relative deviation of 0.1. For `design_var_2`, the actual Uniform distribution depends on the mean value selected by the optimizer for each evaluated design.



