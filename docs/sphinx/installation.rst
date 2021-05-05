.. _installationlabel:

Installation
============

RHEIA is a Python package. The following sections provide information on how to install Python, followed by the installation guide of RHEIA
and the package dependencies for performing deterministic design optimization, robust design optimization and uncertainty quantification.

Installing Python
-----------------

Python can be installed in several ways on your system. If the distribution platform is no constraint,
we recommend installing Python via the `Anaconda Python distribution <https://www.anaconda.com/products/individual>`_, as it includes 
the installation of many common packages in data science (and used in RHEIA), such as NumPy and SciPy.

Installing RHEIA
----------------

RHEIA is available on PyPi, and can be downloaded via the `pip <https://pip.pypa.io/en/stable/>`_ package manager.
The following command installs the most recent version of RHEIA and the package dependencies::

	pip install -i https://test.pypi.org/simple/ rheia-vub
	
Specific from a Jupyter Kernel::

	import sys
	!{sys.executable} -m pip install -i https://test.pypi.org/simple/ rheia-vub
	

Package dependencies
--------------------

The RHEIA features require the installation of several packages.

To evaluate the hydrogen-based energy system models:

- Included in Anaconda:
   - Matplotlib
   - NumPy
   - Pandas 
- Other packages:
   - pvlib
   
To perform uncertainty quantification:

- Included in Anaconda:
   - NumPy
   - SciPy
- Other packages:
   - pyDOE

To perform deterministic design optimization:

- Included in Anaconda:
   - NumPy
- Other packages:
   - pyDOE
   - DEAP

To perform robust design optimization:

- Included in Anaconda:
   - NumPy
   - SciPy
- Other packages:
   - pyDOE
   - DEAP

Import what you need
--------------------

RHEIA allows to import the specific tool you need. To run deterministic or robust design optimization::

	import rheia.OPT.optimization as rheia_opt

To perform uncertainty quantification::

	import rheia.UQ.uncertainty_quantification as rheia_uq

To post-process the results::

    import rheia.POST_PROCESS.lib_post_process as rheia_pp
