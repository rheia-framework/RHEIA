.. _lab:uncertaintyquantification:

Run an uncertainty quantification
=================================

The uncertainty quantification procedure provides the mean and standard deviation of the quantity of interest and the Sobol' indices 
by constructing a Polynomial Chaos Expansion (PCE). More information on PCE is available in :ref:`lab:PCE`.
The mean value and the corresponding uncertainty of the model parameters are characterized in :file:`design_space` and :file:`stochastic_space`, respectively.
More information on characterizing these files is available in :ref:`lab:ssdesignspace` and :ref:`lab:ssstochastic_space`, respectively.  
The model returns the values for the quantity of interest through the function :py:meth:`evaluate()` defined in :py:mod:`case_description`.
More information on this Python wrapper is discussed in :ref:`lab:wrapper`. 


The uncertainty quantification dictionary
-----------------------------------------
To run the uncertainty quantification, first the uncertainty quantification module should be imported::

    import rheia.UQ.uncertainty_quantification as rheia_uq

To characterize the uncertainty quantification, the following dictionary with parameters related to the case and uncertainty quantification should be defined::

    dict_uq = {'case':                  case_name,
               'pol order':             pol_order,
               'objective names':       obj_names,
               'objective of interest': obj_of_interest,
               'results dir':           directory      

               'sampling method':       sampling_method,        #optional, default is 'SOBOL'
               'create only samples':   only_samples_bool,      #optional, default is False
               'draw pdf cdf':          [draw_bool, n_samples], #optional, default is [False]
               'n jobs':                n_jobs,                 #optional, default is 1
               }  

The items of the uncertainty quantification dictionary are described in the following subsections. 
This dictionary is used as the argument for the :py:func`run_uq` function, 
which initiates the uncertainty quantification procedure::

    rheia_uq.run_uq(dict_uq)

Necessary items
^^^^^^^^^^^^^^^

In the following subsections, the necessary items are described.
If one of these items is not provided, the code will return an error.

'case': case_name
~~~~~~~~~~~~~~~~~

The string :py:data:`case_name` corresponds to the name of the case. 
This name should be equal to the name of the folder that comprises the case, which situates in the folder that contains the cases (i.e. :file:`CASES`). 
To illustrate, if the case is defined in :file:`CASES\\CASE_1`, 
the dictionary includes the following item::

		'case': 'CASE_1'


'pol order': pol_order
~~~~~~~~~~~~~~~~~~~~~~

The polynomial order corresponds to the maximum polynomial degree in the PCE.
The polynomial order is characterized by an integer, e.g. for a polynomial order of 2::

	'pol order': 2

'objective names': obj_names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A PCE class-object is constructed for only 1 quantity of interest. However, the statistical moments for several model outputs can be of interest.
To avoid that for each model output, a new set of model evaluations needs to be performed, different model outputs can be stored for each model evaluation.
The names of the different model outputs can be provided in the list :py:data:`'objective_names'`. 
These names are chosen freely by the user, formatted in a string.
If the model returns 3 outputs, the list can be constructed as::

	'objective names': ['output_1', 'output_2', 'output_3']
 
'objective of interest': obj_of_interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Despite that several outputs can be returned for each model evaluation, only one output can be selected as a quantity of interest for the PCE.
The name of this quantity of interest (:py:data:`obj_of_interest`) should be provided. This name should be present in the list of all the objective names.
To illustrate, if the quantity of interest is :py:data:`'output_2'`, then the item in the dictionary is configurated as::

	'objective of interest': 'output_2'

'results dir': directory
~~~~~~~~~~~~~~~~~~~~~~~~

The results directory corresponds to the folder where the results are stored. 
For an illustrative case :py:data:`'CASE_1'`, the UQ results are saved in the folder :file:`RESULTS\\CASE_1\\UQ\\results_1` 
by initiating the following item in the dictionary::

'results dir': results_1

Optional items
^^^^^^^^^^^^^^

The following items are optional items. If one of these items is not provided in the dictionary, 
a default value will be assigned to the key. If none of these are provided, the optional dictionary
items are defined as follows::

               'sampling method':       'SOBOL',
               'create only samples':   False,
               'draw pdf cdf':          [False],
               'n jobs':                1,

'sampling method': sampling_method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the construction of a PCE, a number of model evaluation are required. The samples for model evaluation can be generated
in two different ways: randomly, or through a Sobol' sequence. 
The random generation is called through the string :py:data:`'RANDOM'`, while the Sobol' sequence is initiated through :py:data:`'SOBOL'`.
The default configuration for generating the samples for PCE is through a Sobol' sequence::

	'sampling method': 'SOBOL'
 
'create only samples': only_samples_bool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, the coupling of the system model with the framework is complex. To avoid this coupling, the samples required to determine the statistical moments
can be generated without initiating a model evaluation. Hence, the framework should only generate the samples. To do so,
the Boolean :py:data:`only_samples_bool` can be set to :py:data:`True`::

	'create only samples': True

However, the default configuration sets the value of :py:data:`'create only samples'` to :py:data:`False`::

	'create only samples': False

Additional information on how to create just the samples is present in :ref:`lab:sscreateonlysamples`.

'draw pdf cdf': [draw_bool, n_samples]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the statistical moments, the data for generation the probability density function (pdf) and cumulative distribution function (cdf) can be generated.
This information can be generated by setting the :py:data:`draw_bool` to True and providing the number of samples evaluated on the PCE :py:data:`n_samples`.
To illustrate, to generate pdf and cdf datapoints based on a Monte Carlo evaluation  on the PCE surrogate with 100,000 samples::

    'draw pdf cdf': [True, 1000000]

In the default configuration, the pdf and cdf are not generated::

    'draw pdf cdf': [False]

Reasonable results can be expected on the pdf and cdf when considering 100,000 samples or more. As these samples are evaluated on the PCE,
the computational cost of generating the pdf and cdf is negligible.

'n jobs': n_jobs
~~~~~~~~~~~~~~~~

The number of parallel processes can be defined by the number of available cores on the CPU. 
The default value corresponds to linear processing::

	'n jobs': 1
	
Alternatively, the number of parallel processes can be retreived through the :py:data:`cpu_count` function from the multiprocessing package.
After importing multiprocessing, the item can be defined by::

    'n jobs': int(multiprocessing.cpu_count()/2)

Example of a dictionary for uncertainty quantification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When combining the examples in the previous section, a configurated uncertainty quantification dictionary with the necessary items looks as follows:

.. code-block:: python
   :linenos:

   import rheia.UQ.uncertainty_quantification as rheia_uq

   dict_uq = {'case': 'CASE_1',
              'pol order': 2,
              'objective names': ['output_1', 'output_2', 'output_3'],
              'objective of interest': 'output_2',
              'results dir': 'results_1'      
              }  

   rheia_uq.run_uq(dict_uq)

Alternatively, an uncertainty quantification dictionary which considers random sampling and generates 100,000 PDF and CDF samples on the PCE surrogate:
 
.. code-block:: python
   :linenos:

   import rheia.UQ.uncertainty_quantification as rheia_uq

   dict_uq = {'case': 'CASE_1',
              'pol order': 2,
              'objective names': ['output_1', 'output_2', 'output_3'],
              'objective of interest': 'output_2',
              'results dir': 'results_1'      
              'sampling method': 'RANDOM',
              'draw pdf cdf': [True, 1000000],                
              }  

   rheia_uq.run_uq(dict_uq)

The post-processing of the results is described in :ref:`lab:uqresults`.
	
.. _lab:sscreateonlysamples:

Create samples for an unconnected model
---------------------------------------

When it is burdensome to connect the system model to the framework, the framework provides the option to generate the random samples for uncertainty quantification without a model evaluation.
As the generation of these random samples is based on the characterization of the uncertainties,
a :file:`design_space` file and :file:`stochastic_space` file have to be defined. 
To generate the samples, use (or make a copy of) the :file:`NO_MODEL` folder in :file:`CASES`.
In this folder, a :py:mod:`case_description` module is present, as well as :file:`design_space` and :file:`stochastic_space`.
The :py:mod:`case_description` is empty, as no model evaluations are required.
In :file:`design_space` and :file:`stochastic_space`, the stochastic design space of interest is defined. 
The samples can be generated as follows:

.. code-block:: python
   :linenos:

   import rheia.UQ.uncertainty_quantification as rheia_uq

   dict_uq = {'case':                  'NO_MODEL',
              'pol order':             2,
              'objective names':       ['output_1', 'output_2', 'output_3'],
              'objective of interest': 'output_2',
              'results dir':           'results_1',      
              'create only samples':   True,                
              }  

   rheia_uq.run_uq(dict_uq)

For this example, the samples are written in :file:`RESULTS\\NO_MODEL\\UQ\\results_1\\samples`. Once these samples are evaluated in the model on an external location,
the results can be added to the :file:`RESULTS\\NO_MODEL\\UQ\\results_1\\samples` file. When the results are added for 'output_1', 'output_2', 'output_3', 
the PCE can be constructed for the three quantities of interest. In that case, the value for 'create only samples' is set back to False (i.e. the default value).
To illustrate, for a PCE on 'output_2':

.. code-block:: python
   :linenos:

   import rheia.UQ.uncertainty_quantification as rheia_uq

   dict_uq = {'case':                  'NO_MODEL',
              'pol order':             2,
              'objective names':       ['output_1', 'output_2', 'output_3'],
              'objective of interest': 'output_2',
              'results dir':           'results_1',      
              }  

   rheia_uq.run_uq(dict_uq)

.. warning::
	Make sure that the result directory is equal to the result directory where the updated :file:`samples` file is saved.
	
Post-processing of the results
------------------------------

The results path depends on the case name (e.g. `CASE_1`), the analysis type (UQ)
and the results directory (e.g. `results_1`), i.e. :file:`\\RESULTS\\CASE_1\\UQ\\results_1`.
In this folder, at least 1 file is present: the :file:`samples`  file. This file includes the samples 
and the corresponding deterministic model response, when a system model is connected to the framework (i.e. 'create only samples' set to False).
The second file and third file are named based on the selected maximum polynomial degree and the quantity of interest 
(e.g. :file:`full_pce_order_2_output_2` and :file:`full_pce_order_2_output_2_Sobol_indices` for a polynomial order 2 PCE for the quantity of interest `output_2`).
These files respectively include the PCE information (LOO error, mean and standard deviation) and the Sobol indices (first order and total order).

The Sobol' indices can be represented in a bar chart:

.. code-block:: python
   :linenos:

   import rheia.POST_PROCESS.lib_post_process as rheia_pp

   case = 'case_name'

   pol_order = 1

   my_post_process_uq = rheia_pp.PostProcessUQ(case, pol_order)

   result_dir = 'sample_0'

   objective = 'output_2'

   names, sobol = my_post_process_uq.get_sobol(result_dir, objective)

   plt.barh(names, sobol)
   plt.show()

The LOO-error can be extracted:

.. code-block:: python
   :lineno-start: 17

   loo = my_post_process_uq.get_loo(result_dir, objective)
	
If the data for the Probability Density Function (PDF) and Cumulative Distribution Function (CDF) was generated, both functions can be plotted as follows:

.. code-block:: python
   :lineno-start: 18

   x_pdf, y_pdf = my_post_process_uq.get_pdf(result_dir, objective)

   x_cdf, y_cdf = my_post_process_uq.get_pdf(result_dir, objective)
 


