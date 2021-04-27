

.. _lab:APIref:

API reference
=============
  
Optimization
------------

A function is present to evaluate if the input dictionary is properly
characterized.

.. autosummary::
   :toctree: generated/

   rheia.CASES.determine_stoch_des_space.check_dictionary

The main function that initiates the optimization procedure.

.. autosummary::
   :toctree: generated/
   
   rheia.OPT.optimization.run_opt

In this function, the starting samples are created.

.. autosummary::
   :toctree: generated/
   
   rheia.OPT.optimization.check_existing_results
   rheia.OPT.optimization.scale_samples_to_design_space
   rheia.OPT.optimization.write_starting_samples
   rheia.OPT.optimization.create_starting_samples
      
In addition, the name of the optimization class is loaded.

.. autosummary::
   :toctree: generated/
   
   rheia.OPT.optimization.parse_available_opt
   rheia.OPT.optimization.load_optimizer
   rheia.OPT.genetic_algorithms.return_opt_methods
   rheia.OPT.genetic_algorithms.return_opt_obj   

Finally, an object from the optimization class is instantiated
and the :py:meth:`run_optimizer` is called.
The NSGA-II optimization class includes methods to perform NSGA-II.

.. autosummary::
   :toctree: generated/
   
   rheia.OPT.genetic_algorithms.NSGA2
   rheia.OPT.genetic_algorithms.NSGA2.nsga2_1iter
   rheia.OPT.genetic_algorithms.NSGA2.run_optimizer

The methods to create and evaluate the samples.

.. autosummary::
   :toctree: generated/

   rheia.OPT.genetic_algorithms.NSGA2.define_samples_to_eval
   rheia.OPT.genetic_algorithms.NSGA2.evaluate_samples
   rheia.OPT.genetic_algorithms.NSGA2.assign_fitness_to_population
   rheia.OPT.genetic_algorithms.NSGA2.read_doe
   rheia.OPT.genetic_algorithms.NSGA2.eval_doe
 
The methods to create and update the result files.

.. autosummary::
   :toctree: generated/

   rheia.OPT.genetic_algorithms.NSGA2.init_opt
   rheia.OPT.genetic_algorithms.NSGA2.parse_status
   rheia.OPT.genetic_algorithms.NSGA2.write_status
   rheia.OPT.genetic_algorithms.NSGA2.append_points_to_file

Uncertainty Quantification
--------------------------

The main function that initiates the uncertainty quantification procedure.

.. autosummary::
   :toctree: generated/

   rheia.UQ.uncertainty_quantification.run_uq

This function instantiates an object from :py:class:`Data`. This class includes
methods to acquire the characteristics of the stochastic parameters and to create
the file where the samples are stored

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.Data
   rheia.UQ.pce.Data.create_samples_file
   rheia.UQ.pce.Data.read_stoch_parameters

An object from the class :py:class:`RandomExperiment` is instantiated. This class 
includes a method to determine the number of samples required to construct the PCE.

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.RandomExperiment
   rheia.UQ.pce.RandomExperiment.n_terms

In addition, methods to create the distributions, generate the samples and evaluate 
the samples are present. 

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.RandomExperiment.read_previous_samples
   rheia.UQ.pce.RandomExperiment.create_distributions
   rheia.UQ.pce.RandomExperiment.create_samples
   rheia.UQ.pce.RandomExperiment.create_only_samples
   rheia.UQ.pce.RandomExperiment.evaluate

The PCE class enables to construct a PCE.

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.PCE
   rheia.UQ.pce.PCE.n_to_sum
   rheia.UQ.pce.PCE.multindices
   rheia.UQ.pce.PCE.ols
   rheia.UQ.pce.PCE.calc_a
   rheia.UQ.pce.PCE.run
	
The statistics, Sobol' indices and Leave-One-Out error 
are extracted out of the PCE in the methods below.

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.PCE.get_statistics
   rheia.UQ.pce.PCE.get_psi_sq
   rheia.UQ.pce.PCE.calc_sobol
   rheia.UQ.pce.PCE.calc_loo
	
Finally, the results are printed and stored in corresponding
result files.

.. autosummary::
   :toctree: generated/

   rheia.UQ.pce.PCE.print_res
   rheia.UQ.pce.PCE.draw

To screen the design space (i.e. generate a PCE for a set of design samples),
the following functions allows to retrieve the bounds for the design variables,
to determine the design samples and to generate :file:`design_space` files 
to store the input for the different design samples to be evaluated.

.. autosummary::
   :toctree: generated/

   rheia.UQ.uncertainty_quantification.get_design_variables
   rheia.UQ.uncertainty_quantification.set_design_samples
   rheia.UQ.uncertainty_quantification.write_design_space

Post-processing
---------------

The optimization results are extracted with the methods in :py:class:`PostProcessOpt`.

.. autosummary::
   :toctree: generated/

   rheia.POST_PROCESS.lib_post_process.PostProcessOpt
   rheia.POST_PROCESS.lib_post_process.PostProcessOpt.determine_pop_gen
   rheia.POST_PROCESS.lib_post_process.PostProcessOpt.get_fitness_values
   rheia.POST_PROCESS.lib_post_process.PostProcessOpt.get_population_values
   rheia.POST_PROCESS.lib_post_process.PostProcessOpt.sorted_result_file
   rheia.POST_PROCESS.lib_post_process.PostProcessOpt.get_fitness_population

The uncertainty quantification results are extracted with the methods in :py:class:`PostProcessUQ`.

.. autosummary::
   :toctree: generated/

   rheia.POST_PROCESS.lib_post_process.PostProcessUQ
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.read_distr_file
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.get_sobol
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.get_pdf
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.get_cdf
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.get_loo
   rheia.POST_PROCESS.lib_post_process.PostProcessUQ.get_max_sobol
   
Characterization of design space and stochastic space
-----------------------------------------------------

In :py:class:`StochasticDesignSpace`, methods are present to retrieve the information
on the design space and stochastic space for the specific case.

.. autosummary::
   :toctree: generated/

   rheia.CASES.determine_stoch_des_space.StochasticDesignSpace
   rheia.CASES.determine_stoch_des_space.StochasticDesignSpace.read_design_space
   rheia.CASES.determine_stoch_des_space.StochasticDesignSpace.read_stochastic_space

In addition, a method is present to attach the objectives to the case, as well as a 
method to convert the input sample into a dictionary.

.. autosummary::
   :toctree: generated/

   rheia.CASES.determine_stoch_des_space.StochasticDesignSpace.attach_objectives
   rheia.CASES.determine_stoch_des_space.StochasticDesignSpace.convert_into_dictionary

The cases
---------

The case of interest is loaded in the :py:func:`load_case` function.

.. autosummary::
   :toctree: generated/

   rheia.CASES.determine_stoch_des_space.load_case 

Power-to-fuel
^^^^^^^^^^^^^

A wrapper function is present to evaluate the power-to-fuel model with the samples
generated by the optimization or uncertainty quantification algorithm. In addition,
a function is present to read in the fixed data, required for each model evaluation. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.case_description.evaluate
   rheia.CASES.H2_FUEL.case_description.set_params

The power-to-fuel model contains a class which creates an object that stores information
on the required data.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.h2_fuel.ReadData
   rheia.CASES.H2_FUEL.h2_fuel.ReadData.load_climate
   rheia.CASES.H2_FUEL.h2_fuel.ReadData.load_parameters

In the :py:class:`Evaluation` class, methods are present to quantify the
photovoltaic array power. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.h2_fuel.Evaluation
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.quantify_mpp
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.photovoltaic

In addition, methods are present to characterize the PEM electrolyzer array operation

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.pemel
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.current_to_mh2
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.polyfit_pemel
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.charge_pemel

The :py:meth:`evaluation` method includes the power management strategy.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.evaluation

After applying the power management strategy, the component lifetimes and the 
system costs are quantified, concluded by a method to print the results.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.lifetime
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.cost
   rheia.CASES.H2_FUEL.h2_fuel.Evaluation.print_results

Power-to-power
^^^^^^^^^^^^^^

A wrapper function is present to evaluate the power-to-power model with the samples
generated by the optimization or uncertainty quantification algorithm. In addition,
a function is present to read in the fixed data, required for each model evaluation. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.case_description.evaluate
   rheia.CASES.H2_POWER.case_description.set_params

The power-to-power model contains a class which creates an object that stores information
on the required data.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.ReadData
   rheia.CASES.H2_POWER.h2_power.ReadData.load_climate
   rheia.CASES.H2_POWER.h2_power.ReadData.load_demand
   rheia.CASES.H2_POWER.h2_power.ReadData.load_parameters

In the :py:class:`Evaluation` class, methods are present to generate 
the electricity price profiles and to quantify the
photovoltaic array power. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.Evaluation
   rheia.CASES.H2_POWER.h2_power.Evaluation.elec_profiles
   rheia.CASES.H2_POWER.h2_power.Evaluation.quantify_mpp
   rheia.CASES.H2_POWER.h2_power.Evaluation.photovoltaic
   rheia.CASES.H2_POWER.h2_power.Evaluation.net_power

In addition, methods are present to characterize the PEM electrolyzer array operation

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.Evaluation.pemel
   rheia.CASES.H2_POWER.h2_power.Evaluation.current_to_mh2
   rheia.CASES.H2_POWER.h2_power.Evaluation.polyfit_pemel
   rheia.CASES.H2_POWER.h2_power.Evaluation.charge_pemel

The fuel cell array.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.Evaluation.pemfc
   rheia.CASES.H2_POWER.h2_power.Evaluation.polyfit_pemfc
   rheia.CASES.H2_POWER.h2_power.Evaluation.charge_pemfc

The :py:meth:`evaluation` method includes the power management strategy.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.Evaluation.evaluation

After applying the power management strategy, the component lifetimes,
self sufficiency ratio and system costs are quantified, concluded by a method to print the results.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_POWER.h2_power.Evaluation.lifetime
   rheia.CASES.H2_POWER.h2_power.Evaluation.self_sufficiency_ratio
   rheia.CASES.H2_POWER.h2_power.Evaluation.cost
   rheia.CASES.H2_POWER.h2_power.Evaluation.print_results

Power-to-mobility
^^^^^^^^^^^^^^^^^

A wrapper function is present to evaluate the power-to-mobility model with the samples
generated by the optimization or uncertainty quantification algorithm. In addition,
a function is present to read in the fixed data, required for each model evaluation. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.case_description.evaluate
   rheia.CASES.H2_MOBILITY.case_description.set_params

The power-to-power model contains a class which creates an object that stores information
on the required data.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.ReadData
   rheia.CASES.H2_MOBILITY.h2_mobility.ReadData.load_climate
   rheia.CASES.H2_MOBILITY.h2_mobility.ReadData.load_parameters

In the :py:class:`Evaluation` class, methods are present to generate 
the demand profiles and to quantify the photovoltaic array power. 

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.demand_profiles
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.quantify_mpp
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.photovoltaic

In addition, methods are present to characterize the PEM electrolyzer array operation

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.pemel
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.current_to_mh2
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.polyfit_pemel
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.mh2_to_power

The compressor module.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.compressor
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.polyfit_pemel_compr

The hydrogen tank and dispenser module.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.tank
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.dispenser

The power management strategy module.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.p_for_inst_demand
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.prod_mh2
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.extract_h2_from_tank

The :py:meth:`evaluation` method evaluates the power management strategy.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.evaluation

After applying the power management strategy, the component lifetimes,
system costs and carbon intensity are quantified, concluded by a method to print the results.

.. autosummary::
   :toctree: generated/

   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.lifetime
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.cost
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.lca
   rheia.CASES.H2_MOBILITY.h2_mobility.Evaluation.print_results

Four-bar truss
^^^^^^^^^^^^^^

A wrapper function is present to evaluate the four-bar truss model with the samples
generated by the optimization or uncertainty quantification algorithm.

.. autosummary::
   :toctree: generated/

   rheia.CASES.FOUR_BAR_TRUSS.case_description.evaluate

The four-bar truss model is present in the :py:mod:`four_bar_truss` module.

.. autosummary::
   :toctree: generated/

   rheia.CASES.FOUR_BAR_TRUSS.four_bar_truss.four_bar_truss

EnergyPLAN
^^^^^^^^^^

An evaluate function is present to evaluate the EnergyPLAN Python wrapper with the samples
generated by the optimization or uncertainty quantification algorithm.

.. autosummary::
   :toctree: generated/

   rheia.CASES.ENERGYPLAN.case_description.evaluate

The Python wrapper for the EnergyPLAN model contains functions to create the input text file,
read out the output text file and to run the command which executes the EnergyPLAN executable file.

.. autosummary::
   :toctree: generated/

   rheia.CASES.ENERGYPLAN.run_energyplan.energyplan
   rheia.CASES.ENERGYPLAN.run_energyplan.create_new_input_file
   rheia.CASES.ENERGYPLAN.run_energyplan.read_output_file
