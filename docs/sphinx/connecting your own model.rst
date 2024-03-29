.. _lab:connectingyourownmodel:

Connecting your own model
=========================

The design optimization and uncertainty quantification algorithms are introduced on hydrogen-based energy systems. 
However, RHEIA algorithms can be applied on any system model, thus allowing to evaluate and connect your own model. 
In this section, first, the Python wrapper is described, which connects the system model with the optimization and uncertainty quantification algorithms.
Thereafter, the connection of a Python-based model and a closed-source model is described.

.. _lab:wrapper:

Python wrapper
--------------

The Python wrapper module :py:mod:`case_description` enables to connect the system model to the uncertainty quantification and optimization algorithms.
The module includes two functions: :py:func:`set_params()` and :py:func:`evaluate()`.

In :py:func:`set_params()`, data can be imported before the model evaluations commence. 
In this way, the computational cost of importing data is spent only once,
before the evaluation loop.
To illustrate, importing data from 3 data files (e.g. :file:`datafile_1.csv`, :file:`datafile_2.csv`, :file:`datafile_3.csv`) 
can be defined as follows in the :py:mod:`case_description` module::

    import pandas as pd
	
    def set_params():
        path = os.path.dirname(os.path.abspath(__file__))

        data1 = pd.read_csv(os.path.join(path, 'datafile_1.csv'))['data1']
        data2 = pd.read_csv(os.path.join(path, 'datafile_2.csv'))['data2']
        data3 = pd.read_csv(os.path.join(path, 'datafile_3.csv'))['data3']

        params = [data1, data2, data3]
		
        return params

The model evaluation is performed in the function :py:func:`evaluate()`. This function is called
by the optimization algorithm or uncertainty quantification algorithm, 
which provides the sample as an argument in the form of an enumerate object :py:obj:`x`.
The keys of the dictionary :py:obj:`x[1]` consist of the model parameter names and design variable names, 
as defined in :file:`design_space.csv` and :file:`stochastic_space.csv`.
The corresponding values are the values generated by the optimization or uncertainty quantification algorithm.

In addition to the enumerate object :py:obj:`x`, the :py:func:`evaluate()` function takes :py:data:`params` as an optional argument.
This list contains the list with imported data, as returned by the function :py:func:`set_params()`. 
When :py:func:`set_params()` does not exist, i.e. when no data needs to be imported before the model evaluations commence, 
:py:data:`params` corresponds to an empty list.
In summary, the :py:func:`evaluate()` function can be defined as follows::

    def evaluate(x_in, params = []):
        
        arguments = params + [x_in[1]]

        obj_1, obj_2 = system_model_evaluation(*arguments)
        
        return obj_1, obj_2

Four-bar truss Python model
---------------------------

To illustrate the connection of a Python-based model to RHEIA, a model of a four-bar truss is adopted.
First, the system is briefly illustrated, followed by the model connection and the optimization commands.

Four-bar truss system description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The four-bar truss is presented below:

.. figure:: images/fbt.png
   :scale: 100 %
   :align: center

   The four-bar truss

The aim is to minimize the volume of the truss and to minimize 
the deflection of the outermost joint by controlling the cross-sectional areas of the bars. 
The volume :math:`V` and deflection :math:`d` are defined as:

:math:`V = L (2A_1 + \sqrt{2} A_2 + \sqrt{A_3} + A_4 )`

:math:`d = F L \left( \dfrac{2}{A_1 E_1} + \dfrac{2 \sqrt{2}}{A_2 E_2} - \dfrac{2 \sqrt{2}}{A_3 E_3} + \dfrac{2}{A_4} \right)`

Where :math:`F,~L,~E,~A` are the exerted force, length, Young's modulus and cross-sectional area, respectively. 
The cross-sectional areas are designed respecting the following bounds:

- :math:`A_1 \in [1,3] ~\mathrm{cm}^2`
- :math:`A_2 \in [\sqrt{2},3] ~\mathrm{cm}^2`
- :math:`A_3 \in [\sqrt{2},3] ~\mathrm{cm}^2`
- :math:`A_4 \in [1,3] ~\mathrm{cm}^2`

And the model parameters are defined as follows:

- :math:`L = 200 ~ \mathrm{cm}^2`
- :math:`F \in \mathcal{N}(10,1) ~ \mathrm{kN}`
- :math:`E_1 \in \mathcal{U}(19000,21000) ~ \mathrm{kN}/ \mathrm{cm}^2`
- :math:`E_2 \in \mathcal{U}(19000,21000) ~ \mathrm{kN}/ \mathrm{cm}^2`
- :math:`E_3 \in \mathcal{U}(19000,21000) ~ \mathrm{kN}/ \mathrm{cm}^2`
- :math:`E_4 \in \mathcal{U}(19000,21000) ~ \mathrm{kN}/ \mathrm{cm}^2`

Conclusively, the system model evaluation is coded as follows::

    def four_bar_truss(x_i):

        vol = x_i['L'] * (2. * x_i['A_1'] + 2.**(0.5) *
                          x_i['A_2'] + x_i['A_3']**(0.5) + x_i['A_4'])

        disp = x_i['F'] * x_i['L'] * (2. / (x_i['A_1'] * x_i['E_1']) +
                                      2. * 2**(0.5) / (x_i['A_2'] * x_i['E_2']) -
                                      2. * 2**(0.5) / (x_i['A_3'] * x_i['E_3']) +
                                      2. / (x_i['A_4'] * x_i['E_4']))

        return vol, disp

Where the function argument :py:data:`x` is a dictionary with values for the model parameters, i.e. :math:`L,~F,~E_1,~E_2,~E_3,~E_4`,
and values for the design variables, i.e. :math:`A_1,~A_2,~A_3,~A_4`.
This function is located in the :py:mod:`four_bar_truss` module.

Connecting the case to the framework
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To connect the model to the optimization and uncertainty quantification framework, a specific folder for the model
should be created in the general :file:`CASES` folder. In the :file:`CASES` folder, a reference folder :file:`CASES\\REF` is present, which includes the necessary
files to characterize and connect a new system model. 
Make a copy of the :file:`REF` folder, paste it in the :file:`CASES` folder and rename it, e.g. into :file:`FOUR_BAR_TRUSS`.
Hence, a new case folder is present: :file:`CASES\\FOUR_BAR_TRUSS`.
This folder includes :file:`design_space.csv`, :file:`stochastic_space.csv` and :py:mod:`case_description`.
Finally, include the :py:mod:`four_bar_truss` module in the folder.
The :file:`FOUR_BAR_TRUSS` folder now includes all necessary files and the structure looks as follows::

    CASES
        FOUR_BAR_TRUSS
            design_space.csv
            stochastic_space.csv
            case_despcription.py
            four_bar_truss.py

In :file:`design_space.csv`, the design variables and model parameters are defined.
More information on the characterization of the design space is presented in :ref:`lab:ssdesignspace`.
In the four-bar truss example, 4 design variables (the cross-sectional areas) and 6 model parameters (the force, length and 4 Young's moduli) exist.
The corresponding :file:`design_space.csv` file for the four-bar truss is configured as::

	A_1 var 1 3
	A_2 var 1.414 3
	A_3 var 1.414 3
	A_4 var 1 3
	L   par 200
	F   par 10
	E_1 par 20000
	E_2 par 20000
	E_3 par 20000
	E_4 par 20000

The uncertainty on the stochastic parameters should be defined in :file:`stochastic_space.csv`.
More information on the uncertainty characterization is described in :ref:`lab:ssstochastic_space`.
The exerted force and the Young's moduli are subject to uncertainty.
The corresponding :file:`stochastic_space.csv` file for the four-bar truss file is configured as::

	F   absolute Gaussian 1
	E_1 absolute Uniform  1000
	E_2 absolute Uniform  1000
	E_3 absolute Uniform  1000
	E_4 absolute Uniform  1000

To evaluate the system model in the optimization and uncertainty quantification algorithm, the model should be 
connected to the algorithms. This connection is established in the module :py:mod:`case_description`.
Additional details on this module are provided in :ref:`lab:wrapper`.
At the top of the file, import the module or function that evaluates your model. In the example of the four-bar truss,
this is performed as follows::

	from rheia.CASES.FOUR_BAR_TRUSS.four_bar_truss import four_bar_truss

After the import, the model can be evaluated. This is done in the :py:func:`evaluate()` function.
The dictionary with the input sample names and values :py:data:`x[1]` is passed directly as an argument to the :py:func:`four_bar_truss` function.  
The :py:func:`four_bar_truss` function returns the values for the optimization objectives, i.e. the volume :math:`V` and deflection :math:`d`.
Conclusively, the :py:func:`evaluate()` function is completed as follows::

    def evaluate(x_in, params = []):
        
        vol, disp = four_bar_truss(x_in[1])

        return vol, disp

Run a design optimization
^^^^^^^^^^^^^^^^^^^^^^^^^

Once the characterization and coupling of the case is completed,
the optimization dictionary can be completed to perform the design optimization. 
To illustrate, for a deterministic design optimization:

.. code-block:: python
   :linenos:
    

   import rheia.OPT.optimization as rheia_opt

   dict_opt = {'case':                'FOUR_BAR_TRUSS',
               'objectives':          {'DET': (-1, -1)}, 
               'stop':                ('BUDGET', 9000),
               'population size':     30,
               'results dir':         'run_1',
              }
    
   rheia_opt.run_opt(dict_opt)

In this dictionary, a deterministic design optimization is specified, for which both objectives should be minimized. The computational budget is set at 9000,
which leads to at least 300 generations with a population size of 30. The number of jobs, crossover probability, mutation probability, eta, starting population
and result printing are adopted from the standard setting and are therefore not specified in the dictionary. 
Similarly, the optimization dictionary for robust design optimization on the mean and standard deviation of the displacement can be characterized as follows:

.. code-block:: python
   :linenos:
    

   import rheia.OPT.optimization as rheia_opt

   dict_opt = {'case':                  'FOUR_BAR_TRUSS',
               'objectives':            {'ROB': (-1, -1)}, 
               'stop':                  ('BUDGET', 90000),
               'population size':       30,
               'pol order':             2,
               'objective names':       ['V', 'd'],
               'objective of interest': ['d'],
               'results dir':           'run_1',
               }
    
   rheia_opt.run_opt(dict_opt)

The details on running optimization or uncertainty quantification are provided in 
:ref:`lab:optimization` and :ref:`lab:uncertaintyquantification`, respectively.

EnergyPLAN closed-source model
------------------------------

`EnergyPLAN <https://www.energyplan.eu/>`_ is a software that evaluates national energy system operation. 
It is used by industry, researchers and policy-makers worldwide. 
The software includes the electricity, heating, cooling, industry and transport sector to characterize, among others, the primary energy consumption and CO2 emissions.
Generally, the software is used through a user-friendly user interface, but it can be called through an external command as well.
When calling EnergyPLAN through an external command, the input parameters are provided and the outputs are written in external text files.

In this tutorial, the EnergyPLAN software is connected to the framework. The specific case is based on `exercise 3 <https://www.energyplan.eu/training/exercises/>`_ in the EnergyPLAN training session. 

..
  Note that the aim of this tutorial is to provide a guide on how to connect closed-source software, which is called through an .exe file, 
  for which the input is provided, and the output written, in external files.
  The aim is not to provide novel results , as the considered case is adopted from a simple exercise provided in the EnergyPLAN training session.

The energy system is characterized as follows:

- Electricity demand: :math:`\mathrm{elec\_demand} \in \mathcal{U}(33.3,35.3) ~ \mathrm{TWh/year}`;
- Condensing coal-fired power plant: 9000 MW;
- On-shore wind power: 4000 MW;
- Off-shore wind power: 3000 MW;
- Annual district heating demand: 27.43 TWh (1.59 TWh district heating oil-boilers, 10 TWh small-scale CHP, 15.84 TWh large-scale CHP extraction plants);
- Decentralised natural-gas fired CHP: 1350 MW, with a thermal efficiency of :math:`\mathrm{CHP\_eff\_ht} \in \mathcal{U}(0.45,0.55)` and electrical efficiency of :math:`\mathrm{CHP\_eff\_el} \in \mathcal{U}(0.36,0.46)`;
- Large-scale coal-fired CHP: 2000 MW, with a thermal efficiency of 50% and electrical efficiency of 41%;
- Heat Pump of 300 MWe with a :math:`\mathrm{COP} \in \mathcal{U}(2.5,3.5)`;
- Individual house Fuel demand for heating:14.42 TWh (0.01 coal, 4.2 oil, 5.66 natural gas and 4.55 biomass);
- Industrial fuel demand: 53.66 TWh (3.37 coal, 26.92 oil, 18.19 natural gas and 5.18 biomass);
- Industrial district heating production: 1.73 TWh;
- Industrial district electricity production: 2.41 TWh;
- Transportation fuel demand: 13.25 TWh Jet Petrol, 27.50 TWh Diesel, 23.35 TWh Petrol and 2.55 TWh hydrogen;
- 900 MWe electrolyzers at :math:`\mathrm{eff\_H2} \in \mathcal{U}(0.6,0.7)` efficiency for hydrogen production for transport.

Connecting the model to the framework
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, a specific folder for the model should be created in the :file:`CASES` folder, e.g. :file:`ENERGYPLAN`.
The necessary files are :file:`design_space.csv`, the :py:mod:`case_description` module( and :file:`stochastic_space.csv` for uncertainty quantification and robust design optimization).
In addition, a module to call the EnergyPLAN software :py:mod:`run_energyplan` 
and a .csv file to provide the input that represent the current case (:file:`case.csv`) are included.
This results in the following structure::

    CASES
        ENERGYPLAN
            design_space.csv
            stochastic_space.csv
            case_despcription.py
            run_energyplan.py
            case.csv

The :file:`design_space.csv` file includes the mean values for the stochastic model parameters::

	elec_demand par 34.3
	COP         par 3
	eff_H2      par 0.65
	CHP_eff_el  par 0.41
	CHP_eff_th  par 0.5

More information on the characterization of the design space is presented in :ref:`lab:ssdesignspace`.
The stochastic space is defined in :file:`stochastic_space.csv`::

	elec_demand absolute Uniform 1
	COP         absolute Uniform 0.5
	eff_H2      absolute Uniform 0.05
	CHP_eff_el  absolute Uniform 0.05
	CHP_eff_th  absolute Uniform 0.05

More information on the uncertainty characterization is described in :ref:`lab:ssstochastic_space`.

In the :py:mod:`run_energyplan` module, the :py:func:`energyplan` Python wrapper evaluates the EnergyPLAN model and returns the model outputs of interest.
In this function, first, the index and sample are splitted from the :py:data:`x_tup` argument. 
Then, the path of the :file:`energyPLAN.exe` executable and the :file:`case.csv` file with input parameters are provided::

    x = x_tup[1]
    index = x_tup[0]
    path = os.path.dirname(os.path.abspath(__file__))

    ep_path = r'C:\energyPLAN\energyPLAN.exe'
    input_file = os.path.join(path, 'case.csv')

A new input file is created for the energyPLAN model, 
where the initial values of the model parameters are overwritten by the new values, provided by the input sample :py:data:`x`::

    new_input_file = '%s_%i.csv' % (input_file[:-4], index)

    create_new_input_file(input_file, new_input_file, x)

This new input file, with updated values on the model parameters for each evaluation, has a name that ends with the index of the input sample, 
i.e. the sample position in the array of samples created for uncertainty quantification.
In this way, a unique input file is generated for each input sample, which ensures that during the parallelization of the model evaluations over the available CPUs,
the input file that corresponds to the input sample is evaluated. Following a similar logic, the name of the file with results is defined as::

    result_file = os.path.join(path, 'result_%i.csv' % index)

Once the EnergyPLAN input file is characterized, the command to run EnergyPLAN is called::

    cm = [ep_path, '-i', new_input_file, '-ascii', result_file]
    sp.call(cm)

Followed by reading the values for the quantities of interest::

    co2, fuel = read_output_file(result_file)

We refer to the :py:mod:`run_energyplan` module for additional details on the :py:func:`create_new_input_file` and :py:func:`read_output_file` functions. 

In the :py:mod:`case_description` module, the function to run the case in EnergyPLAN is imported from the :py:func:`run_energyplan()` module at the top of the script::

    from rheia.CASES.ENERGYPLAN.run_energyplan import energyplan

The :py:func:`run_energyplan()` function is evaluated in :py:func:`evaluate()`, where the enumerate object :py:data:`x_in` is provided as an argument. For this case, no fixed parameters are provided as an argument::

    def evaluate(x_in, params=[]):

        co2, fuel = energyplan(x_in)
        
        return co2, fuel

The enumerate object :py:data:`x_in` contains the sample to be evaluated and the index of the sample in the list of samples.

Run uncertainty quantification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the characterization complete for uncertainty quantification, the algorithm can be initiated with:

.. code-block:: python
   :linenos:
    

   import rheia.UQ.uncertainty_quantification as rheia_uq
   import multiprocessing as mp

	
   dict_uq = {'case':                  'ENERGYPLAN',
              'n jobs':                int(mp.cpu_count() / 2),
              'pol order':             1,
              'objective names':       ['co2', 'fuel'],
              'objective of interest': 'fuel',
              'results dir':           'run_1'      
              }  

   if __name__ == '__main__':
       rheia_uq.run_uq(dict_uq)
	
..
  The results illustrate a LOO-error 0.005 for both the primary energy supply and CO2-emission.
  For illustration purposes, the Sobol' indices for primary energy supply and CO2-emission are shown below.

  .. figure:: images/cyom_sobol_fuel.png
   :scale: 100 %
   :align: center

   The Sobol' indices for the primary energy supply.

  .. figure:: images/cyom_sobol_co2.png
   :scale: 100 %
   :align: center

   The Sobol' indices for the CO2-emission.

The details on running optimization or uncertainty quantification are provided in 
:ref:`lab:optimization` and :ref:`lab:uncertaintyquantification`, respectively.
