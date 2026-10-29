.. _lab:stochasticdesignspace:

Define the design and stochastic spaces
=======================================

In this section, the characterization of the model parameters and design variables is described.
The design variables are variables that are controllable by the designer within a certain range. 
This range should be provided, to shape the search space for the optimizer. 
A model parameter corresponds to a parameter that usually cannot be controlled by the designer (e.g. the cost of a photovoltaic panel), 
or the decision on this parameter value is fixed (e.g. a fixed amount of photovoltaic panels on the roof). 
Such a parameter can be considered deterministic or uncertain.
The design variables and model parameters are characterized in two files: :file:`design_space.csv` and :file:`stochastic_space.csv`.
These files are present in the case folder for the specific case considered (e.g. :file:`CASES\\H2_FUEL` for the `H2_FUEL` case).
In :file:`design_space.csv`, the deterministic values for the model parameters and the range for the design variables are provided.
In :file:`stochastic_space.csv`, the uncertainty is allocated to the specific model parameters and design variables.
When a deterministic design optimization is performed, only the :file:`design_space.csv` file is required.
In the other cases, i.e. uncertainty quantification and robust design optimization, both files are required.
When the uncertainty on a parameter is known as two percentiles rather than as a mean and a deviation,
both files can be generated from those estimates instead of being filled in by hand (see :ref:`lab:sscharacterizeuncertainty`).

.. _lab:ssdesignspace:

The design_space.csv file
-------------------------

In the :file:`design_space.csv` file, the design variables and the model parameters which need a quantification in your model are defined. 
When performing uncertainty quantification, the :file:`design_space.csv` file consists only of model parameters.
In the case of deterministic design optimization or robust design optimization, the :file:`design_space.csv` file requires design variables. 
Additionally, if some model parameters require a quantification outside the model, 
the :file:`design_space.csv` file includes both :ref:`lab:ssdesignvariables` and :ref:`lab:ssmodelparameters`.

.. _lab:ssdesignvariables:

Characterizing the design variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
To define a design variable, the set-up in the :file:`design_space.csv` file is as follows::

	name,type,value,upper_bound

where:

- name: name of the variable;
- type: a design variable is indicated with `var`;
- value: lower bound value for the design variable;
- upper_bound: upper bound value for the design variable. 

An example for a configured design variable `des_var_1` with a range between 10 and 50 is::

    des_var_1,var,10,50


.. _lab:ssmodelparameters:

Characterizing the model parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this file :file:`design_space.csv`, the deterministic value (or mean value when the parameter is considered uncertain) is provided.
The configuration of a model parameter is similar to the configuration of a design variable::

    name,type,value,upper_bound

where:

- name: name of the variable;
- type: a parameter is indicated with `par`;
- value: deterministic value (or mean value when the parameter is stochastic).
- upper_bound: left empty for model parameters.

An example of a configured model parameter `par_1` with a mean value of 0.03 is::

	par_1,par,0.03,

.. _lab:ssexampleds:

Example of design_space.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^
An example of a configured :file:`design_space.csv` file, which consists of 3 model parameters (par_1, par_2 and par_3) and 2 design variables (design_var_1 and design_var_2), is presented::

	name,type,value,upper_bound
	design_var_1,var,1,3
	design_var_2,var,10,100
	par_1,par,4,
	par_2,par,2.5,
	par_3,par,175,

.. _lab:ssstochastic_space:

The stochastic_space.csv file
-----------------------------

The uncertainty on the design variables and model parameters can be allocated through the file :file:`stochastic_space.csv`. 
This file is required when performing robust design optimization and uncertainty quantification, where several parameters are subjected to uncertainty. 
For every design variable and model parameter defined in :file:`design_space.csv`, an uncertainty can be defined.

Characterizing the uncertainties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Defining the uncertainty can be done as follows::

	name,relation,distribution,deviation

where:

	- name: name of the parameter or variable, equal to the name of the parameter or variable in :file:`design_space.csv`;
	- relation: absolute or relative uncertainty to the mean, defined with `absolute` or `relative`, respectively;
	- distribution: The distribution of the uncertainty;
	- deviation: uncertainty on the mean.

The following distributions are available:

- Uniform
- Gaussian
- Lognormal

The meaning of deviation at the end of the line depends on the distribution. When a Uniform distribution is considered,
the deviation refers to the absolute (or relative) difference between the upper bound of the Uniform distribution and the mean: for :math:`\mathcal{U}(a,b)`, :math:`deviation = (b-a)/2`). 
When a Gaussian distribution is considered, the value corresponds to the standard deviation: :math:`\mathcal{N}(mean,deviation)`.
When a Lognormal distribution is considered, the mean and deviation are the physical mean and standard deviation of the lognormal variable.
The generated samples in :file:`samples.csv` and the values passed to the model are lognormal physical values.
For PCE construction, RHEIA transforms these values to the corresponding latent standard Gaussian variable and uses Hermite polynomials.
Keep always in mind that the mean value of the design variables is the deterministic value assigned by the optimizer in each iteration. In the case of model fixed parameter the mean value corresponds to the fixed value as it is assigned in the model definition.

An example of a configured uncertain parameter `par_2`, characterized by a Uniform distribution with a :math:`\pm 1` deviation from the mean value::

	par_2,absolute,Uniform,1

Note that it is not required to allocate an uncertainty to every design variable and model parameter defined in :file:`design_space.csv`.
In other words, when a parameter (or variable) is defined in :file:`design_space.csv`, but not in :file:`stochastic_space.csv`, the parameter (or variable) is considered deterministic. 
Moreover, it is not necessary to keep the same order of appearance of parameters and variables in :file:`design_space.csv` :file:`stochastic_space.csv` files.

Example of stochastic_space.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In summary, a :file:`stochastic_space.csv` file corresponding to the illustrative :file:`design_space.csv` example file in :ref:`lab:ssexampleds` looks like this::

	name,relation,distribution,deviation
	par_1,relative,Gaussian,0.5
	par_2,absolute,Uniform,1
	design_var_2,relative,Uniform,0.1

Where the model parameter `par_3` and design variable `design_var_1` are considered deterministic,
`par_1` is characterized by a Gaussian distribution with a
relative standard deviation of 0.5 (i.e. :math:`\mathcal{N}(4,2)`),
`par_2` is characterized by a Uniform distribution with an
absolute deviation of 1 (i.e. :math:`\mathcal{U}(1.5,3.5)`) and
`design_var_2` is characterized by a Uniform distribution with a
relative deviation of 0.1. For `design_var_2`, the actual Uniform distribution depends on the mean value selected by the optimizer for each evaluated design.

.. _lab:sscharacterizeuncertainty:

Characterizing uncertainty from two percentiles
-----------------------------------------------

The two files described above ask for a mean and a deviation.
Uncertainty, however, is often available as two percentiles instead,
e.g. a declared value at the 50th percentile and a design value at the 99th percentile from a standard or a specification.
Two percentiles pin down a two-parameter distribution exactly, with no other data needed,
so such estimates can be converted directly into the mean and deviation that :file:`design_space.csv` and :file:`stochastic_space.csv` expect.
For a Lognormal input, the fit follows Cook, J.D. (2010), *Determining distribution parameters from quantiles*.

The module :py:mod:`rheia.UQ.characterize_uncertainty` performs this conversion for the three distributions supported by RHEIA's PCE
(Uniform, Gaussian and Lognormal), and writes the result straight into a case's :file:`design_space.csv` and :file:`stochastic_space.csv`.

.. _lab:sscharacterizationfile:

The uncertainty_characterization.csv file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The percentile estimates are collected in a :file:`uncertainty_characterization.csv` file,
placed in the case folder next to :file:`design_space.csv` and :file:`stochastic_space.csv`.
A template with the expected header is available in :file:`CASES\\NO_MODEL\\uncertainty_characterization.csv`, ready to be copied into your own case folder.
Every uncertain parameter is defined on one line::

	name,distribution,x_low,p_low,x_high,p_high

where:

- name: name of the parameter;
- distribution: the distribution to fit, i.e. `Uniform`, `Gaussian` or `Lognormal`;
- x_low: the value estimated at the lower percentile;
- p_low: that lower percentile, as a fraction between 0 and 1;
- x_high: the value estimated at the upper percentile;
- p_high: that upper percentile, as a fraction between 0 and 1.

Both percentiles should satisfy :math:`0 < p_{low} < p_{high} < 1`, with :math:`x_{low} < x_{high}`.
A Lognormal fit additionally requires both values to be strictly positive.

An example of a configured file, for a material's thermal conductivity declared as 500 W/mK at the 50th percentile
and 3000 W/mK at the 99th percentile, together with a capital cost estimated between 400 and 550 at the 5th and 95th percentile::

	name,distribution,x_low,p_low,x_high,p_high
	th_cond,Lognormal,500,0.5,3000,0.99
	capex,Uniform,400,0.05,550,0.95

Running the characterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the case folder contains the characterization file, the conversion is performed with :py:func:`characterize_from_csv`::

	import rheia.UQ.characterize_uncertainty as rheia_char

	rheia_char.characterize_from_csv('MY_CASE')

For every row, the fitted mean is written into :file:`design_space.csv` as a model parameter (`par`),
and the fitted deviation into :file:`stochastic_space.csv` with an `absolute` relation,
since that is what two absolute-valued percentiles determine.
Both files are created when they do not exist yet, so a case folder containing nothing but the characterization file is enough to get started.
In that situation, a :py:mod:`case_description` module is copied from the :file:`NO_MODEL` template as well (see :ref:`lab:sscreateonlysamples`),
so the case works immediately with :py:data:`'create only samples': True`,
and :py:func:`evaluate()`/:py:func:`set_params()` can be filled in later to connect a real model.

When the files already exist, the parameters are merged into them:
a parameter that is not characterized yet is appended, while a parameter that is already present is replaced by the updated estimate.
Rows that the characterization file does not mention are left untouched, so an existing case can be extended or refreshed without disturbing the rest of it.
Which parameters were added and which were replaced is printed, so an update is never silent::

	Characterized 2 parameter(s) in case 'MY_CASE' (...\CASES\MY_CASE):
	  added:   capex
	  updated: th_cond
	  written to design_space.csv and stochastic_space.csv

Passing :py:data:`overwrite=False` raises on a parameter that already exists, instead of replacing it.

Iterating on the percentiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Picking good percentiles is usually iterative.
Setting :py:data:`write=False` fits the distributions without touching :file:`design_space.csv` or :file:`stochastic_space.csv`,
so the same call can be repeated freely while tweaking the estimates.
Combined with :py:data:`plot=True` and :py:data:`show=True`, the fitted probability density function is rendered with the two elicited percentiles marked,
as a sanity check on the fit::

	# iterate: look at the fit, nothing written yet
	rheia_char.characterize_from_csv('MY_CASE', write=False, plot=True, show=True)

	# happy with the fit: write it for real
	rheia_char.characterize_from_csv('MY_CASE')

Once :py:data:`write=True`, the plot of each parameter is saved as a :file:`<name>_fit.png` file in the case folder.
A complete, editable example is provided in :file:`examples\\characterize_uncertainty.py`.

Finally, the percentile estimates can also be passed directly in Python, rather than through a file,
with :py:func:`characterize_parameters`::

	rheia_char.characterize_parameters(
	    'MY_CASE',
	    [{'name': 'th_cond', 'distribution': 'Lognormal',
	      'x_low': 500, 'p_low': 0.50, 'x_high': 3000, 'p_high': 0.99}])



