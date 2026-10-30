.. _lab:whatsnew:

This page contains a description of the latest updates of RHEIA.

What's new
==========

v2.3.1.
-------

- The columns of :file:`samples.csv` are identified by the names on its
  header line. The header was written but never read back, so every column
  was taken from a fixed position: a samples file filled in externally with
  its quantities of interest in another order than :py:data:`'objective
  names'` was read from the wrong column, without any error. The statistics
  and Sobol' indices then described a different model output than the one
  requested. See :ref:`how the columns are identified
  <lab:sssamplescolumns>`.
- When the header does not name every stochastic parameter and quantity of
  interest exactly once, the columns are read in the order they appear, as
  before, and a warning reports what was expected and what was found. Files
  written by earlier versions are therefore read exactly as they were.
- Adding samples to a file whose columns are in a different order than the
  run dictionary is refused. New rows are written in the order of the run
  dictionary, so appending would otherwise leave two column orders in one
  file.

v2.3.0.
-------

This release makes runs reproducible, fixes three defects that could affect
users directly, and adds continuous integration.

Reproducibility
^^^^^^^^^^^^^^^

- A new optional :py:data:`'seed'` item in the optimization and uncertainty
  quantification dictionaries makes a run reproducible. It covers the
  starting population, the crossover and mutation operators, and the samples
  drawn for the Polynomial Chaos Expansion. See :ref:`lab:ssoptseed` and
  :ref:`lab:ssuqseed`.
- :py:func:`set_design_samples` accepts a :py:data:`seed` argument for the
  same reason.
- Omitting the item leaves every sampler behaving exactly as it did before,
  so existing results are unchanged.

Bug fixes
^^^^^^^^^

- A run dictionary can be reused. :py:func:`run_opt` attached the evaluation
  function to the caller's own dictionary, which was then rejected on the
  next call, so a script that ran two optimizations from one dictionary
  failed on the second.
- A starting population containing a design sample that sits exactly on a
  design variable bound is loaded correctly. Such samples are produced by
  the crossover and mutation operators, so restarting from a previous
  population could fail.
- The :file:`examples/run_uq_h2_fuel.py` example carried the objective names
  of the four-bar truss case, so it reported hydrogen production under the
  label :code:`displacement`.
- The warning that DEAP emits when it overwrites its :py:data:`Fitness` and
  :py:data:`Individual` classes, which happens on every optimization after
  the first in a session, is no longer shown.

Documentation
^^^^^^^^^^^^^

- The optimization dictionary no longer documents an :py:data:`'algo'` item.
  That item was rejected as an unknown key, and the alternative optimizers it
  referred to are not available; NSGA-II is the optimizer RHEIA provides.

Development
^^^^^^^^^^^

- A continuous integration workflow runs the test suite on Python 3.11, 3.12
  and 3.13 for Linux and Windows, together with a style check.

v2.2.0.
-------

This release adds a step to characterize uncertain model parameters from
percentile estimates, so the design space and stochastic space can be
generated from elicited values rather than filled in by hand.

Uncertainty characterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- A new module :py:mod:`rheia.UQ.characterize_uncertainty` converts two
  percentile estimates of a model parameter into the mean and deviation that
  :file:`design_space.csv` and :file:`stochastic_space.csv` expect. Two
  percentiles pin down a two-parameter distribution exactly; the Lognormal fit
  follows Cook, J.D. (2010), *Determining distribution parameters from
  quantiles*.
- The estimates are provided in a new :file:`uncertainty_characterization.csv`
  file in the case folder, with one row per uncertain parameter
  (:code:`name,distribution,x_low,p_low,x_high,p_high`). A template is shipped
  in :file:`CASES\\NO_MODEL`.
- :py:func:`characterize_from_csv` reads that file and writes the fitted
  parameters into the case. Both space files are created when absent, so a case
  folder holding only the characterization file is enough to get started; a
  :py:mod:`case_description` module is then copied from the :file:`NO_MODEL`
  template, so the case runs immediately with
  :py:data:`'create only samples': True`.
- Parameters that are not characterized yet are appended, while parameters that
  are already present are replaced by the updated estimate. Rows that the
  characterization file does not mention are left untouched. Which parameters
  were added and which were replaced is printed, so an update is never silent.
- :py:data:`write=False` fits without touching disk, and :py:data:`plot=True`
  renders the fitted probability density function with the two elicited
  percentiles marked, to support iterating on the estimates before committing
  them.
- The percentile estimates can also be passed directly in Python with
  :py:func:`characterize_parameters`, instead of through a file.
- The supported distributions are the same three the PCE supports: Uniform,
  Gaussian and Lognormal.
- An :file:`examples/characterize_uncertainty.py` example was added.

v2.1.0.
-------

This release extends uncertainty quantification, improves optimization
post-processing and makes sample generation for external models easier to use.

Uncertainty quantification
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Lognormal stochastic inputs are now supported in :file:`stochastic_space.csv`
  with :code:`Lognormal` as distribution type.
- Lognormal samples written to :file:`samples.csv` and passed to the model
  remain physical lognormal values. Internally, RHEIA transforms them to the
  corresponding latent Gaussian variables for Hermite PCE construction.
- Existing :file:`samples.csv` files with lognormal physical samples are
  transformed in the same way before fitting the PCE.
- Model evaluations are appended to :file:`samples.csv` one by one and flushed
  immediately, so completed evaluations are preserved if a run stops early.
- UQ now validates that the number of outputs returned by the model matches
  the number of entries in :py:data:`'objective names'`.
- A :file:`NO_MODEL` create-only workflow example was added to generate PCE
  training samples without a connected model.

Optimization and post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- :py:meth:`PostProcessOpt.get_hypervolume` calculates the hypervolume of the
  Pareto front for every generation.
- Hypervolume post-processing supports mixed minimization and maximization
  objectives through objective weights.
- The deterministic optimization tutorial now includes a hypervolume
  convergence plot.

Testing
^^^^^^^

- Added workflow coverage for create-only PCE training samples with
  :file:`NO_MODEL`.
- Added an end-to-end lognormal equivalence test comparing a lognormal input
  model against an equivalent latent-Gaussian model.

v2.0.0.
-------

This release modernizes the package and input formats and adds sparse polynomial chaos expansion support.

Uncertainty quantification
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Sparse PCE is now available through the uncertainty quantification dictionary with :py:data:`'uq method': 'sparse'`.
- Sparse PCE training size is configured with :py:data:`'n samples'`.
- Full PCE remains the default through :py:data:`'uq method': 'full'`.
- UQ result files now distinguish full and sparse PCE results. Sparse files include the number of training samples in the filename, e.g. :file:`sparse_pce_order_2_lcoh_n_samples_80.txt`.
- The post-processing class :py:class:`PostProcessUQ` accepts :py:data:`method` and :py:data:`samples` to read either full or sparse PCE results.
- The PCE implementation uses SciPy's maintained quasi-Monte Carlo tools for Sobol' sequences.

Input files
^^^^^^^^^^^

- :file:`design_space.csv` now has a required header: :code:`name,type,value,upper_bound`.
- :file:`stochastic_space.csv` now has a required header: :code:`name,relation,distribution,deviation`.
- Parameter rows in :file:`design_space.csv` keep :code:`upper_bound` empty; variable rows require both :code:`value` and :code:`upper_bound`.
- All packaged case files and tutorial design-space files have been migrated to the headered schema.
- Reading and writing of design-space and stochastic-space CSV files now uses pandas and validates required columns explicitly.

Optimization and post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Latin Hypercube Sampling now uses SciPy's maintained sampler.
- Optimization post-processing uses pandas for sorted CSV output and validates generation consistency in :file:`fitness.csv` and :file:`population.csv`.
- UQ post-processing uses pandas for PDF, CDF and Sobol' CSV files and reports clearer errors for malformed PCE summaries.

Packaging
^^^^^^^^^

- Package metadata and dependencies moved to :file:`pyproject.toml`.
- Legacy :file:`setup.py` and :file:`manifest.in` were removed.
- RHEIA now requires Python 3.11 or newer.
- The documentation build installs the local checkout instead of the published PyPI package.
- A :file:`.gitignore` was added for common Python, build, IDE and generated result files.

v1.1.11.
--------

This was the previous stable version of RHEIA. 

v1.1.5.
-------

The input files for the cases (:file:`design_space.csv` and :file:`stochastic_space.csv`)
and the output files for optimization (:file:`fitness.csv`, :file:`population.csv`)
and uncertainty quantification (:file:`samples.csv`, :file:`full_pce_order_x_y_Sobol_indices.csv`) are converted from extensionless files into CSV files.

v1.0.0.
-------

This is the official release of RHEIA! It corresponds to the package published in the Journal of Open Source Software.
