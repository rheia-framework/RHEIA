"""
The :py:mod:`post_process` module contains classes to present the
optimization and uncertainty quantification results.
"""

import os
import numpy as np


class PostProcessOpt():
    """
    The PostProcessOpt class provides methods to retrieve the fitness and
    population for the generations provided during the optimization.

    Parameters
    ----------
    case : str
        Name of the case to be evaluated.
    eval_type : str
        The type of evaluation performed. Equals 'DET' for deterministic
        optimization and 'ROB' for robust optimization.

    """

    def __init__(self, case, eval_type):
        path_file = os.path.dirname(os.path.abspath(__file__))
        path_start = os.path.abspath(os.path.join(path_file, os.pardir))
        path = os.path.join(path_start,
                                 'RESULTS',
                                 case
                                 )
        self.n_pop = 0
        self.x_lines = []
        self.y_lines = []
        self.result_path = os.path.join(path,
                                        eval_type,
                                        )
        self.fitness_file = ''
        self.population_file = ''
        self.n_gen = 0

    def determine_pop_gen(self):
        """

        Determines the number of design samples in the population
        and the number of generations performed.

        """

        self.n_gen = 0
        with open(self.fitness_file, 'r') as file:
            self.y_lines = file.readlines()
            for string in self.y_lines:
                if string == '- \n':
                    # count the number of generations
                    self.n_gen += 1

            # determine the population size
            self.n_pop = int(len(self.y_lines) / self.n_gen - 1)

    def get_fitness_values(self, gen):
        """

        Returns the fitness values for the population
        generated in the specified generation.

        Parameters
        ----------
        gen : int
            The generation of interest.

        Returns
        -------
        fit_val : ndarray
            The fitness values for the population of interest.

        """

        fit_val = np.zeros((len(self.y_lines[0].split()), self.n_pop))
        for index, line in enumerate(
                self.y_lines[(gen - 1) * self.n_pop - 1 + gen:
                             gen * self.n_pop - 1 + gen]):

            # add fitness values for the desired generation
            fit_val[:, index] = [float(i) for i in line.split()]

        return fit_val

    def get_population_values(self, gen):
        """
        Returns the fitness values for the population
        generated in the specified generation.

        Parameters
        ----------
        gen : int
            The generation of interest.

        Returns
        -------
        pop_val : ndarray
            The population of interest.

        """

        with open(self.population_file, 'r') as file:
            self.x_lines = file.readlines()

        pop_val = np.zeros((len(self.x_lines[0].split()), self.n_pop))
        for index, line in enumerate(
                self.x_lines[(gen - 1) * self.n_pop - 1 + gen:
                             gen * self.n_pop - 1 + gen]):

            # add population values for the desired generation
            pop_val[:, index] = [float(i) for i in line.split()]

        return pop_val

    def sorted_result_file(self, y_val, x_val):
        """

        Generates the files that include the sorted population
        and fitness files. The population and fitness are sorted
        based on the first objective.

        Parameters
        ----------
        y_val : ndarray
            The fitness values.
        x_val : ndarray
            The set of design samples.

        """

        # write sorted fitness values in corresponding file
        with open(self.fitness_file + "_final_sorted", 'w') as file:
            for sample in y_val:
                for value in sample:
                    file.write('%f, ' % value)
                file.write('\n')

        # write sorted population values in corresponding file
        with open(self.population_file + "_final_sorted", 'w') as file:
            for sample in x_val:
                for value in sample:
                    file.write('%f, ' % value)
                file.write('\n')

    def get_fitness_population(self, result_dir, gen=0):
        """

        Returns the population and corresponding fitness values
        for the generation of interest.

        Parameters
        ----------
        result_dir : str
            The directory were the results are stored.
        gen : int, optional
            The generation of interest. The default is 0,
            i.e. the final generation.

        Returns
        -------
        y_gen : ndarray
            The fitness values.
        x_gen : ndarray
            The set of design samples.

        """

        self.fitness_file = os.path.join(self.result_path,
                                         result_dir,
                                         'fitness',
                                         )

        self.population_file = os.path.join(self.result_path,
                                            result_dir,
                                            'population',
                                            )

        # set the population size and number of generations
        self.determine_pop_gen()

        # get fitness values for desired generation
        y_unsorted = self.get_fitness_values(gen)

        # get population values for desired generation
        x_unsorted = self.get_population_values(gen)

        a, b = x_unsorted.shape
        c, d = y_unsorted.shape

        # get indices when sorting the y_unsorted array,
        # sorting is based on first objective
        indices = np.argsort(y_unsorted[0])

        x_gen = np.zeros((a, b))
        y_gen = np.zeros((c, d))
        for j, k in enumerate(indices):
            for index, y_in in enumerate(y_gen):
                y_gen[index][j] = y_unsorted[index][k]
            for index, x_in in enumerate(x_gen):
                x_gen[index][j] = x_unsorted[index][k]

        # print sorted results in separate files
        self.sorted_result_file(y_gen.transpose(), x_gen.transpose())

        return y_gen, x_gen


class PostProcessUQ():
    """
    The PostProcessUQ class provides methods to retrieve the LOO error,
    plot the Sobol indices, PDF and CDF.

    Parameters
    ----------
    case : str
        Name of the case to be evaluated.
    pol_order : int
        The polynomial order.

    """

    def __init__(self, case, pol_order):
        path_file = os.path.dirname(os.path.abspath(__file__))
        path_start = os.path.abspath(os.path.join(path_file, os.pardir))
        path = os.path.join(path_start,
                                 'RESULTS',
                                 case
                                 )
        self.result_path = os.path.join(path,
                                        'UQ',
                                        )

        self.pol_order = pol_order

    def read_distr_file(self, distr_file):
        """

        Reads the file with information on the
        cumulative density function or probability
        density function.

        Parameters
        ----------
        distr_file : str
            The name of the distribution file.

        Returns
        -------
        x_val : ndarray
            The values from the PDF or CDF on the
            quantity of interest.
        y_val : ndarray
            The probability density (for the PDF)
            or cumulative probability (for the CDF).

        """

        # read the file with info on the CDF or PDF
        with open(distr_file, 'r') as file:
            lines = file.readlines()
            x_val = np.ones(len(lines) - 1)
            y_val = np.ones(len(lines) - 1)
            for index, line in enumerate(lines[1:]):
                tmp = line.split()
                x_val[index] = float(tmp[0])
                y_val[index] = float(tmp[1])

        return x_val, y_val

    def get_sobol(self, result_dir, objective):
        """

        Retrieves the information on the Sobol' indices from
        the corresponding file in the result directory.

        Parameters
        ----------
        result_dir : str
            The result directory.
        objective : str
            The name of the quantity of interest.

        Returns
        -------
        names : list
            The names of the stochastic parameters.
        sobol : list
            The total' order Sobol' indices.

        """

        sobol_file = os.path.join(self.result_path,
                                  '%s' % result_dir,
                                  'full_pce_order_%i_%s_Sobol_indices' % (
                                      self.pol_order, objective)
                                  )

        # retrieve the parameter names and corresponding Sobol indices
        res_tmp = []
        with open(sobol_file, 'r') as file:
            for line in file.readlines()[1:]:
                res_tmp.append([i for i in line.split()])
            names = [row[0] for row in res_tmp]
            sobol = [float(row[2]) for row in res_tmp]

        return names, sobol

    def get_pdf(self, result_dir, objective):
        """

        Retrieves the points that define the probability density function.

        Parameters
        ----------
        result_dir : str
            The result directory.
        objective : str
            The name of the quantity of interest.

        Returns
        -------
        x_val : ndarray
            The values from the PDF on the
            quantity of interest.
        y_val : ndarray
            The probability density.

        """

        pdf_file = os.path.join(self.result_path,
                                '%s' % result_dir,
                                'data_pdf_%s' % objective
                                )

        # get the x and y values for the PDF
        x_val, y_val = self.read_distr_file(pdf_file)

        return x_val, y_val

    def get_cdf(self, result_dir, objective):
        """

        Retrieves the points that define the cumulative density function.

        Parameters
        ----------
        result_dir : str
            The result directory.
        objective : str
            The name of the quantity of interest.

        Returns
        -------
        x : ndarray
            The values from the CDF on the
            quantity of interest.
        y : ndarray
            The cumulative probability.

        """

        cdf_file = os.path.join(self.result_path,
                                '%s' % result_dir,
                                'data_cdf_%s' % objective
                                )

        # get the x and y values for the CDF
        x_val, y_val = self.read_distr_file(cdf_file)

        return x_val, y_val

    def get_loo(self, result_dir, objective):
        """

        Reads the Leave-One-Out error from the corresponding
        file in the result directory.

        Parameters
        ----------
        result_dir : str
            The result directory.
        objective : str
            The name of the quantity of interest.

        Returns
        -------
        loo : float
            The Leave-One-Out error.

        """

        loo_file = os.path.join(self.result_path,
                                '%s' % (result_dir),
                                'full_pce_order_%i_%s' % (
                                    self.pol_order, objective)
                                )

        # retrieve the LOO error
        with open(loo_file, 'r') as file:
            line = file.readlines()[0]
            loo = float(line.split()[1])

        return loo

    def get_max_sobol(self, result_dirs, objective, threshold=0.05):
        """

        This method gathers the Sobol' indices for each stochastic parameter
        for each sample. The highest Sobol' index for each stochastic
        parameter is compared with the threshold value. If the highest
        Sobol' index is higher than the threshold, the name of the
        stochastic parameter is printed under 'significant Sobol indices'.
        If not, it is printed under 'negligible Sobol indices'.

        Parameters
        ----------
        result_dir : list
            The result directories.
        objective : str
            The name of the quantity of interest.
        threshold : float, optional
            The threshold that determines if a Sobol' index
            is considered significant. The default is 0.05.

        """

        # store the dictionary with parameter names and Sobol indices in a
        # list for each result directory evaluated
        n_samples = len(result_dirs)
        res_dict = [{}] * n_samples
        for index, result_dir in enumerate(result_dirs):
            names, sobol = self.get_sobol(result_dir, objective)
            res_dict[index] = dict(zip(names, sobol))

        # get the highest Sobol index for each parameter over the different
        # dictionaries
        max_dict = dict()
        for name in names:
            sobol_res = np.zeros(n_samples)
            for j, dic in enumerate(res_dict):
                sobol_res[j] = dic[name]
            max_dict[name] = max(sobol_res)

        print('significant Sobol indices:')
        for k in names:
            if max_dict[k] >= threshold:
                print('%s: %4f' % (k, max_dict[k]))

        print('\nnegligible Sobol indices:')
        for k in names:
            if max_dict[k] < threshold:
                print('%s: %4f' % (k, max_dict[k]))
