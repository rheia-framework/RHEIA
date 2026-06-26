.. _lab:whatsnew:

This page contains a description of the latest updates of RHEIA.

What's new
==========

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
