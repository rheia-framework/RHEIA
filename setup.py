import setuptools




NAME = "rheia"
DESCRIPTION = 'Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems'
URL = 'https://github.com/rheia-framework/RHEIA'
EMAIL = 'rheia.framework@gmail.com'
AUTHOR = 'Diederik Coppitters, Panagiotis Tsirikoglou, Ward De Paepe, Konstantinos Kyprianidis, Anestis Kalfas, Francesco Contino'
VERSION = '1.1.11'

with open("README.md", "r", encoding="utf-8") as fh:
    LONG_DESCRIPTION = fh.read()


setuptools.setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description = LONG_DESCRIPTION,
      url=URL,
      author=AUTHOR,
      author_email=EMAIL,
      packages= setuptools.find_packages(),
      classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
      ],       
      install_requires=[
      'pyDOE>=0.3.8',
      'deap==1.3.1',
      'numpy>=1.24.1',
      'scipy>=1.10.0',
      'sobolsequence>=0.2.1',
      'pandas>=1.5.3',
      'matplotlib>=3.2.2',
      'pvlib>=0.9.4',
      'h5py>=3.8.0'
      ],
      python_requires = ">=3.6",
      include_package_data=True,
      zip_safe=False)
