import setuptools




NAME = "rheia"
DESCRIPTION = 'Robust design optimization of renewable Hydrogen and dErIved energy cArrier systems'
URL = 'https://github.com/rheia-framework/RHEIA'
EMAIL = 'rheia.framework@gmail.com'
AUTHOR = 'Diederik Coppitters, Panagiotis Tsirikoglou, Ward De Paepe, Konstantinos Kyprianidis, Anestis Kalfas, Francesco Contino'
VERSION = '1.0.0'

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
      'pyDOE',
      'deap',
      'numpy',
      'scipy',
      'sobolsequence'
      ],
      python_requires = ">=3.6",
      include_package_data=True,
      zip_safe=False)