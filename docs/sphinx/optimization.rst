.. _lab:optimization:

Run an optimization task
========================

The deterministic design optimization procedure optimizes model outputs of interest by searching a finite design space as it is constructed by selected model input parameters, also known as design variables. 
The robust design optimization works under the same principle. The fundamental change is focused on the probabilistic treatment of the model input parameters (i.e. definition and propagation of uncertainties), which dictate the optimization of mean and minimization of the standard deviation of the considered model outpus.
The multi-objective optimization algorithm used in RHEIA is Nondominating Sorting Genetic Algorithm (NSGA-II). More information on NSGA-II is available in :ref:`lab:ssnsga2`.
The design variables and model parameters are characterized in the :file:`design_space.csv` file.
As per robust design optimization, the uncertainty on the respective design variables and model parameters is characterized in the :file:`stochastic_space.csv` file.
More information on characterizing these files is available in :ref:`lab:ssdesignspace` and :ref:`lab:ssstochastic_space`, respectively. 
The system model evaluations are coupled with the optimization algorithm in :py:mod:`case_description.`.
More information on this Python wrapper is discussed in :ref:`lab:wrapper`. 
 

.. _lab:ssrundetopt:

Run deterministic design optimization
---------------------------------------

To run a deterministic design optimization, first the optimization module should be imported::

    import rheia.OPT.optimization as rheia_opt

To characterize and configure a deterministic design optimization, the following Python dictionary should be completed::

    dict_opt = {'case':                case_name,
                'objectives':          {opt_type: weights}, 
                'population size':     n_pop,
                'stop':                ('BUDGET', comp_budget),
                'results dir':         directory,

                'x0':                  (pop_type, pop_method), #optional, default is ('AUTO', 'LHS')
                'cx prob':             c_prob,                 #optional, default is 0.9
                'mut prob':            mut_prob,               #optional, default is 0.1
                'eta':                 eta,                    #optional, default is 0.2
                'n jobs':              n_jobs,                 #optional, default is 1 
                }

This dictionary is used as the argument for the :py:func:`run_opt()` function, 
located in the :py:mod:`optimization` module, which starts the optimization procedure::

    rheia_opt.run_opt(dict_opt)

.. warning::
   When parallel processing is considered (i.e. :py:data:`dict_opt['n jobs']` > 1), 
   calling the :py:func:`run_opt()` function should be protected as the main 
   entry point using :py:data:`if __name__ == '__main__'`.

The dictionary includes necessary items and optional items. The items are clarified in the following sections.

Necessary items
^^^^^^^^^^^^^^^

In the following subsections, the necessary items are described.
If one of these items is not provided, the code will return an error.

'case': case_name
~~~~~~~~~~~~~~~~~

The string :py:data:`case_name` corresponds to the name of the case. 
This name should be equal to the name of the folder that enclosures all the case files, which situates in the folder that contains the cases (i.e. :file:`CASES`). 
To illustrate, if the optimization case is defined in :file:`CASES\\CASE_1`, 
the dictionary includes the following item::

		'case': 'CASE_1'

'objectives': {opt_type: (weights)} 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the item with :py:data:`'objectives'` key, the optimization type and the weigths for the objectives are specified. 
Two optimization types are available: deterministic design optimization (:py:data:`'DET'`) and robust design optimization (:py:data:`'ROB'`).
The weights are defined in a tuple and determine if the objective is either maximized or minimized.
When minimization of an objective is desired, the weight corresponds to -1. 
Instead, when maximization is desired, the weight corresponds to 1. 
For deterministic design optimization (:py:data:`'DET'`), the order of the weights corresponds to the order of the model outputs
returned by the method :py:meth:`evaluate()` (see :ref:`lab:wrapper`).  
For instance, for 2 objectives which should be minimized simultaneously in a deterministic design optimization, the dictionary item reads::

	'objectives': {'DET': (-1, -1)}

Alternatively, maximizing the first objective and minimizing the second and third objective corresponds to::

	'objectives': {'DET': (1, -1, -1)}
	
In the robust design optimization approach, the mean and standard deviation for each quantity of interest is optimized.
For each quantity of interest, the weight for the mean and standard deviation should be provided.
Hence, the weights with even index correspond to the mean, while the weights with odd index correspond to the standard deviation.
To illustrate, when the mean should be maximized and the standard deviation minimized for two quantities of interest, the dictionary item reads::

	'objectives': {'ROB': (1, -1, 1, -1)}

Instead, when only one quantity of interest is desired, for which both the mean and standard deviation should be minimized, the item reads::

	'objectives': {'ROB': (-1, -1)}
	
Note that for robust design optimization, the number of weights should be equal to two times the number of quantities of interest (i.e. the mean and standard deviation for each
quantity of interest is an objective). Therefore, make sure that the number of quantities of interest defined (see :ref:`lab:secobjofint`) matches the number of weights defined.

'population size': n_pop
~~~~~~~~~~~~~~~~~~~~~~~~~~

The population size corresponds to the number of samples contained in a single population. 
After each evaluation of the entire population, the optimizer generates a new population with an equal number of samples.
This iterative process continues until the predefined computational budget is complied with. 
Hence, with a computational budget of 1440 model evaluations, 
a population size of 20 will lead to at least 72 generations for deterministic design optimization::

	'population size': 20
	
Note that when the population number and computational budget do not result in an integer for the number of generations, 
the number of generations is rounded up to the nearest integer.  
Additional details on defining the value for the population size is illustrated in :ref:`lab:choosepop`. 

'stop': ('BUDGET', comp_budget)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stopping criterion for the optimization is defined by the computational budget, i.e. the number of model evaluations. 
This is a common engineering stopping criterion, which is defined based on the time available
to perform the optimization. To illustrate, when the system model takes 10 seconds to evaluate and 4 cores are available for parallel processing, 
the computational budget for a deterministic design optimization procedure of 1 hour is equal to 1440. 
The allocation of this computational budget through the integer :py:data:`comp_budget` is illustrated below::

	'stop': ('BUDGET', 1440)

'results dir': directory
~~~~~~~~~~~~~~~~~~~~~~~~

The result directory corresponds to the folder where the results are stored. 
For an illustrative deterministic design optimization (:py:data:`'DET'`) of a case (:py:data:`'CASE_1'`), the results are stored in the folder :file:`RESULTS\\CASE_1\\DET\\results_1` 
by initiating the following key-value pair in the dictionary::

'results dir': 'results_1'

.. warning::
	If previous results exist in the results directory, the optimization procedure continues from the last, previously generated, population. 
	Hence, any specified characterization of the initial population in the optimization dictionary is ignored. However, the computational budget is renewed. 

.. _lab:optitemsdet:

Optional items
^^^^^^^^^^^^^^

In addition to the necessary items, optional items can be added to the dictionary. 
If one of these items is not provided in the dictionary, a typical value will be assigned to the key. 
The default configuration for these optional items is::

                'x0':                  ('AUTO', 'LHS'), 
                'cx prob':             0.9,
                'mut prob':            0.1,
                'eta':                 0.2,
                'n jobs':              1, 

.. _lab:ssx0:

'x0': (pop_type, pop_method) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Information can be provided to characterize the starting population. If no information is available on the starting population, 
the population can be generated automatically by defining the string :py:data:`pop_type` with :py:data:`'AUTO'`. 
When :py:data:`'AUTO'` is selected, there are two ways of generating the population automatically: 
randomly (:py:data:`pop_method` = :py:data:`'RANDOM'`) or through Latin Hypercube Sampling (:py:data:`pop_method` = :py:data:`'LHS'`). 
The default configuration for this item is the generation of the first population through LHS::

	'x0': ('AUTO', 'LHS')

Alternatively, when information on the starting population is available, the :py:data:`pop_type` should be defined by :py:data:`'CUSTOM'`. 
In that case, the starting population should be provided in a separate file,
located in the case folder. The name of the file corresponds to the string that defines :py:data:`pop_method`. 
To illustrate for :py:data:`'CASE_1'`, with a starting population saved in :file:`CASES\\CASE_1\\x0_start.csv`, the item is defined as::

	'x0': ('CUSTOM', 'x0_start.csv')

This CSV file should contain a number of samples equal to the population size. 
Each sample is characterized by a number of values equal to the number of design variables, delimited by a white space.
Each value should situate between the lower bound and upper bound of the corresponding design variable, 
in the order of appearance of the design variables in the :file:`design_space.csv` file.

Example: 

The following design variables are defined in :file:`design_space.csv`::

	var_1,var,1,3
	var_2,var,0.4,0.9
	var_3,var,12,15

Then, for a population size of 5, a suitable characterization of the starting population file is::

	1.43,0.78,13.9,
	2.97,0.44,12.1,
	1.12,0.64,14.2,
	2.31,0.51,14.5,
	2.05,0.88,13.6,
	
For instance, the optimized population retrieved from a previous run can be copied from the :file:`population_final_sorted.csv` file, 
and pasted in the :file:`x0_start.csv` file.

'cx prob': c_prob
~~~~~~~~~~~~~~~~~

The probability of crossover at the mating of two parent samples.
The default crossover probability is equal to 0.9::

	'cx prob': 0.9
	
More information on setting the crossover probability is illustrated in :ref:`lab:choosepop`. 

'mut prob': mut_prob
~~~~~~~~~~~~~~~~~~~~

The probability of mutation, i.e. the probability of values in the design samples being flipped.
The default value on the mutation probability corresponds to::

	'mut prob': 0.1

More information on setting the mutation probability is illustrated in :ref:`lab:choosepop`. 

'eta': eta
~~~~~~~~~~

The crowding degree of the crossover, which determines the resemblance of the children to their parents. 
The default crowding degree is::

    'eta': 0.2

'n jobs': n_jobs
~~~~~~~~~~~~~~~~

The number of parallel processes can be defined by the number of available cores on the Central Processing Unit. 
The default value corresponds to linear processing::

	'n jobs': 1
	
Alternatively, the number of parallel processes can be retreived through the :py:data:`cpu_count` function from the multiprocessing package.
After importing multiprocessing, the item can be defined by::

    'n jobs': int(multiprocessing.cpu_count()/2)

Example of a dictionary for deterministic design optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When combining the examples in the previous section, a fully-defined optimization dictionary with the necessary items looks as follows:

.. code-block:: python
   :linenos:
    

   import rheia.OPT.optimization as rheia_opt

   dict_opt = {'case':                'CASE_1',
               'objectives':          {'DET': (-1, -1)}, 
               'population size':     20,
               'stop':                ('BUDGET', 1440),
               'results dir':         'results_1',
               }

   rheia_opt.run_opt(dict_opt)

In the example below, parallel processing is considered, the optimization starts from a predefined population, defined in :py:data:`'x0_start.csv'`, 
and the crossover probability is decreased to 0.85:

.. code-block:: python
   :linenos:

   import rheia.OPT.optimization as rheia_opt
   import multiprocessing as mp

   dict_opt = {'case':                'CASE_1',
               'objectives':          {'DET': (-1, -1)}, 
               'population size':     20,
               'stop':                ('BUDGET', 1440),
               'results dir':         'results_1',
               'x0':                  ('CUSTOM', 'x0_start.csv'), 
               'cx prob':             0.85,
               'n jobs':              int(mp.cpu_count() / 2),
               }

   if __name__ == '__main__':
       rheia_opt.run_opt(dict_opt)

.. _lab:runrdo:

Run a robust design optimization
--------------------------------

For robust design optimization, like for deterministic design optimization, first the optimization module should be imported::

    import rheia.OPT.optimization as rheia_opt

To characterize the robust design optimization, the following dictionary with parameters related to the case, optimization 
and uncertainty quantification should be provided::

    dict_opt = {'case':                  case_name,
                'objectives':            {opt_type: weights}, 
                'population size':       n_pop,
                'stop':                  ('BUDGET', comp_budget),
                'results dir':           directory,
                'pol order':             pol_order,
                'objective names':       obj_names,
                'objective of interest': obj_of_interest,

                'x0':                    (pop_type, pop_method), #optional, default is ('AUTO', 'LHS')
                'cx prob':               c_prob,                 #optional, default is 0.9
                'mut prob':              mut_prob,               #optional, default is 0.1
                'eta':                   eta,                    #optional, default is 0.2
                'n jobs':                n_jobs,                 #optional, default is 1 
                'sampling method':       sampling_method         #optional, default is 'SOBOL'
                }

This dictionary is used as the argument for the :py:func:`run_opt()` function, which starts the optimization procedure::

    rheia_opt.run_opt(dict_opt)

The necessary and optional items in the dictionary for deterministic design optimization are also present in the dictionary for robust design optimization.
These items are described in :ref:`lab:ssrundetopt`.
The additional necessary and optional items for robust design optimization are described in the following subsections. 

Necessary items
^^^^^^^^^^^^^^^

In the following subsections, the additional necessary items for robust design optimization are described.
If one of these items is not provided, the code will return an error.


'pol order': pol_order
~~~~~~~~~~~~~~~~~~~~~~

The polynomial order corresponds to the maximum polynomial degree in the PCE trunctation scheme.
The polynomial order is characterized by an integer, e.g. for a polynomial order of 2::

	'pol order': 2
	
Determining the appropriate polynomial order is strongly case-specific. A method to determine the order is presented in the next section :ref:`lab:detpolorder`.

'objective names': [obj_names]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model might return several outputs (i.e. for multi-objective optimization).
The names of the different model outputs can be provided in the list :py:data:`objective_names`. 
These names are chosen freely by the user, formatted in a string.
If the model returns 3 outputs, the list can be constructed as::

	'objective names': ['output_1', 'output_2', 'output_3']
 

.. _lab:secobjofint:

'objective of interest': obj_of_interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Despite that several outputs can be returned for each model evaluation, not all outputs might be of interest for the robust design optimization.
The quantities of interest should be provided in the list :py:data:`obj_of_interest`. These names should be present in the list of all the objective names.
To illustrate, for a robust design optimization with the mean and standard deviation of :py:data:`'output_2'` and :py:data:`'output_3'` as objectives, 
the item in the dictionary is configured as::

	'objective of interest': ['output_2','output_3']

Instead, if a robust design optimization is desired with :py:data:`'output_3'` as quantity of interest::

	'objective of interest': ['output_3']

Optional items
^^^^^^^^^^^^^^

When running robust design optimization, only one additional optional item exists, in addition to the 
optional items presented in the deterministic design optimization section (:ref:`lab:optitemsdet`).
The item is described below.

'sampling method': sampling_method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the construction of a PCE, a number of model evaluation are required (see :ref:`lab:pce`). These samples can be generated
in two different ways: randomly, or through a Sobol' sequence. 
The random generation is called through the string :py:data:`'RANDOM'`, while the Sobol' sequence is initiated through :py:data:`'SOBOL'`.
The default configuration for generating the samples for PCE is through a Sobol' sequence::

	'sampling method': 'SOBOL'

Example of a dictionary for robust design optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When combining the examples in the previous section, a configurated optimization dictionary with only necessary items for robust design optimization looks as follows:

.. code-block:: python
   :linenos:

   import rheia.OPT.optimization as rheia_opt

   dict_opt = {'case':                  'CASE_1',
               'objectives':            {'ROB': (-1,-1,-1,-1)}, 
               'population size':       20,
               'stop':                  ('BUDGET', 1440),
               'results dir':           'results_1',
               'pol order':             2,
               'objective names':       ['output_1', 'output_2', 'output_3'],
               'objective of interest': ['output_2','output_3'],
               }

   rheia_opt.run_opt(dict_opt)

An additional example, where parallel processing is considered, the mutation probability is decreased to 0.05 and the sampling method is random:

.. code-block:: python
   :linenos:

   import rheia.OPT.optimization as rheia_opt
   import multiprocessing as mp

   dict_opt = {'case':                  'CASE_1',
               'objectives':            {'ROB': (-1,-1,-1,-1)}, 
               'population size':       20,
               'stop':                  ('BUDGET', 1440),
               'results dir':           'results_1',
               'pol order':             2,
               'objective names':       ['output_1', 'output_2', 'output_3'],
               'objective of interest': ['output_2','output_3'],
               'mut prob':              0.05,
               'sampling method':       'RANDOM',
               'n jobs':                int(mp.cpu_count()/2), 
               }

   if __name__ == '__main__':
       rheia_opt.run_opt(dict_opt)

The post-processing of the results is described in :ref:`lab:optimizationresults`.

.. _lab:detpolorder:

Screening of the design space
-----------------------------

Considering the current truncation scheme, the polynomial order and the number of stochastic parameters define the number of model evaluations 
required to construct the PCE (see :ref:`labpce`). In robust design optimization, a PCE is constructed for each design sample evaluated during the optimization.
Hence, the polynomial order should be sufficient over the entire design space. In addition, only the stochastic parameters which have a significant
impact on the standard deviation on the quantity of interest. To determine the polynomial order and the significant stochastic parameters, a screening of
the design space is performed as follows:

.. code-block:: python
   :linenos:

   import rheia.UQ.uncertainty_quantification as rheia_uq
   import multiprocessing as mp

   case = 'case_name'    
   var_dict = rheia_uq.get_design_variables(case)

   n_samples = 5
   X = rheia_uq.set_design_samples(var_dict, n_samples)
    
   for iteration, x in enumerate(X):

       rheia_uq.write_design_space(case, iteration, var_dict, x)

       dict_uq = {'case':                  case,
                  'pol order':             1,
                  'objective names':       ['obj_1','obj_2'],
                  'objective of interest': 'obj_1',
                  'results dir':           'sample_%i' %iteration      
                  }  

       rheia_uq.run_uq(dict_uq, design_space = 'design_space_%i.csv' %iteration)

After providing the name of the case, a dictionary with the design variable names, lower bounds and upper bounds can be defined
via the :py:func:`get_design_variables` function.
From this dictionary, the design samples can be constructed through LHS via :py:func:`set_design_samples`. 
Then, for each design sample in the array :py:data:`X`, a :file:`design_space.csv` file is constructed through the function :py:func:`write_design_space()`. 
For each :file:`design_space.csv` file, the PCE is constructed through the characterization of the uncertainty quantification dictionary. 
For more information on the characterization of this dictionary, we refer to :ref:`lab:uncertaintyquantification`.
The uncertainty quantification dictionary and the specific :file:`design_space.csv` file is then provided to the :py:func:`run_uq` function.
This results in a PCE for each design sample, with a corresponding Leave-One-Out (LOO) error. That LOO error is stored in the :file:`RESULTS` folder.
Considering the specific dictionary determined above, the results for the different design samples are stored in :file:`\\RESULTS\\case_name\\UQ`::

    RESULTS 
        case_name
            UQ
                sample_0
                sample_1
                sample_2
                sample_3
                sample_4
	
Where in each folder, the LOO error is stored in :file:`full_PCE_order_1_obj_1`.

Determine the polynomial order
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The maximum polynomial degree for the multivariate polynomials needs to be determined up front and its value should ensure accurate
statistical moments on the quantity of interest in the considered stochastic space. An indication on the accuracy of the PCE is
the Leave-One-Out (LOO) error. If the error is below a certain threshold, the PCE achieves an acceptable accuracy. This threshold is a user-defined constant. 
To ensure accurate statistical moments during the robust design optimization procedure, the polynomial order should be sufficient 
over the entire design space. In other words, for each design sample, the polynomial order should be sufficient to construct an accuracte PCE.
Latin Hypercube Sampling is used to construct a set of design samples, which provides a representation of the design space. If the worst-case LOO 
among the corresponding PCEs is still below a certain threshold, the corresponding polynomial order can be considered sufficient to be used during
the robust design optimization procedure.

The worst-case LOO error (i.e. the highest LOO error over the diffferent design samples) can be determined as follows:

.. code-block:: python
   :linenos:

   import rheia.POST_PROCESS.post_process as rheia_pp

   case = 'case_name'

   pol_order = 1
   my_post_process_uq = rheia_pp.PostProcessUQ(case, pol_order)

   n_samples = 5
   result_dirs = ['sample_%i' %i for i in range(n_samples)]

   objective = 'obj_1'

   loo = [0]*n_samples
   for index, result_dir in enumerate(result_dirs):
       loo[index] = my_post_process_uq.get_loo(result_dir, objective)

   print(max(loo))

Where the :py:meth:`get_loo()` method returns the LOO error for every sample.
Based on the worst-case LOO error, the maximum polynomial degree of the PCE for the robust design optimization can be evaluated.

Reducing the stochastic dimension
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The contribution of each parameter uncertainty to the variance of the quantity of interest is different. 
When the stochastic parameters that contribute little to the output variance are considered deterministic,
the computational efficiency can be improved dramatically, with a negligible loss in accuracy on the statistical moments. 
To make sure that the Sobol' index for a specific parameter is negligible over the entire design space, 
the Sobol' indices from the screening are adopted. 
The highest Sobol' index found for each stochastic parameter over the set of design samples
determines the Sobol' index on which the decision is made in this conservative approach.
The stochastic parameters with negligible effect are printed through the following commands, 
where a threshold for the Sobol' index is set at 1/number of uncertain parameters (in this example, 10 uncertain parameters):

.. code-block:: python
   :linenos:

   import rheia.POST_PROCESS.post_process as rheia_pp

   case = 'case_name'

   pol_order = 1
   my_post_process_uq = rheia_pp.PostProcessUQ(case, pol_order)

   n_samples = 5
   result_dirs = ['sample_%i' %i for i in range(n_samples)]

   objective = 'obj_1'

   my_post_process_uq.get_max_sobol(result_dirs, objective, threshold=1./10.)	

.. warning::
	As the accuracy of this method depends mainly on the number of design samples considered, the results are only indicative.
	Therefore, the stochastic parameters with negligible Sobol' index are not removed automatically. It is suggested to evaluate the feasibility of
	this result, based on the knowledge of the user on the considered system model.
	
Post-processing of the results
------------------------------

An illustrative path directs towards the result files from optimization, 
for which the path depends on the case name (e.g. :py:data:`'CASE_1'`), the analysis type (:py:data:`'DET'` or :py:data:`'ROB'`)
and the results directory (e.g. :py:data:`'results_1'`), is defined as follows: :file:`\\RESULTS\\CASE_1\\DET\\results_1`.
In this folder, 3 folder are present: :file:`STATUS.txt`, :file:`fitness.csv` and :file:`population.csv`.
The :file:`STATUS.txt` file consists of two columns: ITER and EVALS. In ITER, the finished generation number is saved, while the corresponding number in EVALS
provides the actual computational budget spent after completing that generation.
The :file:`population.csv` and :file:`fitness.csv` files contain the design samples and results, respectively. 
This information is stored for every design sample in every generation. 
The design sample on line :math:`j` in :file:`population.csv` corresponds to the fitness 
on line :math:`j` in :file:`fitness.csv`.
Plotting the results can be performed as follows:

.. code-block:: python
   :linenos:

   import rheia.POST_PROCESS.post_process as rheia_pp
   import matplotlib.pyplot as plt

   case = 'case_name'

   eval_type = 'DET'

   my_opt_plot = rheia_pp.PostProcessOpt(case, eval_type)

   result_dir = 'run_1'

   y, x = my_opt_plot.get_fitness_population(result_dir)

   plt.plot(y[0], y[1], '-o')
   plt.xlabel('obj_1')
   plt.ylabel('obj_2')
   plt.show()

   for x_in in x:
       plt.plot(y[0], x_in, '-o')
   plt.legend(['x_1', 'x_2'])
   plt.xlabel('obj_1')
   plt.ylabel('x')
   plt.show()

The method :py:meth:`get_fitness_population()` returns, for the last available generation, the fitness values and the population.
Alternatively, a number of generations can be plotted on the same graph by defining the optional argument :py:data:`gen`. 
This enables to evaluate the convergence of the result. To illustrate, plotting 
generation 5, 15 and 25 can be done as follows:

.. code-block:: python
   :lineno-start: 25
	
   for i in [5,15,25]:
       y,x = my_opt_plot.get_fitness_population(result_dir, gen = i)
       plt.plot(y[0], y[1])
   plt.xlabel('obj_1')
   plt.ylabel('obj_2')
   plt.show()

When calling the :py:meth:`get_fitness_population()` method, the design samples and fitness values are sorted based on the first objective and saved in :file:`population_final_sorted.csv` 
and :file:`fitness_final_sorted.csv`, respectively, in the results directory.
