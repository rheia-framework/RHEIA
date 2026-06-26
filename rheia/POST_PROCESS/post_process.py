"""
The :py:mod:`post_process` module contains classes to present the
optimization and uncertainty quantification results.
"""

import os
import numpy as np
import pandas as pd


class PostProcessOpt:
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
        self._fitness_generations = []
        self._population_generations = []

    @staticmethod
    def _read_generation_file(filename):
        """
        Read a generation-separated CSV file into one array per generation.
        """
        generations = []
        current = []

        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                if line == '-':
                    if current:
                        generations.append(np.asarray(current, dtype=float))
                        current = []
                    continue

                current.append([float(value) for value in line.split(',')
                                if value.strip() != ''])

        if current:
            generations.append(np.asarray(current, dtype=float))

        return generations

    def determine_pop_gen(self):
        """

        Determines the number of design samples in the population
        and the number of generations performed.

        """

        self._fitness_generations = self._read_generation_file(self.fitness_file)
        self._population_generations = self._read_generation_file(
            self.population_file)

        self.n_gen = len(self._fitness_generations)
        if self.n_gen == 0:
            raise ValueError("No generations found in %s." % self.fitness_file)
        if self.n_gen != len(self._population_generations):
            raise ValueError(
                "Fitness and population files contain a different number "
                "of generations.")

        self.n_pop = len(self._fitness_generations[0])
        for generation, population in zip(self._fitness_generations,
                                          self._population_generations):
            if len(generation) != self.n_pop or len(population) != self.n_pop:
                raise ValueError("Population size changes between generations.")

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

        if gen <= 0:
            gen = self.n_gen
        if gen > self.n_gen:
            raise ValueError("Generation %i does not exist." % gen)

        return self._fitness_generations[gen - 1].transpose()

    def get_population_values(self, gen):
        """
        Returns the design variables for the population
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

        if gen <= 0:
            gen = self.n_gen
        if gen > self.n_gen:
            raise ValueError("Generation %i does not exist." % gen)

        return self._population_generations[gen - 1].transpose()

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

        pd.DataFrame(y_val).to_csv(
            self.fitness_file[:-4] + "_final_sorted.csv",
            header=False, index=False, lineterminator='\n')
        pd.DataFrame(x_val).to_csv(
            self.population_file[:-4] + "_final_sorted.csv",
            header=False, index=False, lineterminator='\n')

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
                                         'fitness.csv',
                                         )

        self.population_file = os.path.join(self.result_path,
                                            result_dir,
                                            'population.csv',
                                            )

        # set the population size and number of generations
        self.determine_pop_gen()

        # get fitness values for desired generation
        y_unsorted = self.get_fitness_values(gen)

        # get population values for desired generation
        x_unsorted = self.get_population_values(gen)

        # get indices when sorting the y_unsorted array,
        # sorting is based on first objective
        indices = np.argsort(y_unsorted[0])

        y_gen = y_unsorted[:, indices]
        x_gen = x_unsorted[:, indices]

        # print sorted results in separate files
        self.sorted_result_file(y_gen.transpose(), x_gen.transpose())

        return y_gen, x_gen

    @staticmethod
    def _nondominated(points):
        """
        Return the nondominated points for a minimization problem.
        """
        points = np.unique(np.asarray(points, dtype=float), axis=0)
        keep = np.ones(len(points), dtype=bool)

        for i, point in enumerate(points):
            if not keep[i]:
                continue
            dominated = np.all(points <= point, axis=1) & np.any(
                points < point, axis=1)
            if np.any(dominated):
                keep[i] = False

        return points[keep]

    @staticmethod
    def _hypervolume_minimization(points, reference_point):
        """
        Calculate the dominated hypervolume for minimization objectives.
        """
        points = np.asarray(points, dtype=float)
        reference_point = np.asarray(reference_point, dtype=float)

        if len(points) == 0:
            return 0.

        if np.any(points > reference_point):
            raise ValueError(
                "The reference point should be worse than all Pareto points.")

        if points.shape[1] == 1:
            return reference_point[0] - np.min(points[:, 0])

        bounds = np.unique(np.concatenate((points[:, 0],
                                           [reference_point[0]])))
        bounds.sort()

        volume = 0.
        for i in range(len(bounds) - 1):
            width = bounds[i + 1] - bounds[i]
            if width <= 0.:
                continue

            active = points[points[:, 0] <= bounds[i], 1:]
            if len(active) == 0:
                continue

            volume += width * PostProcessOpt._hypervolume_minimization(
                active, reference_point[1:])

        return volume

    def get_hypervolume(self, result_dir, reference_point,
                        objective_weights=None):
        """
        Calculate the hypervolume of the Pareto front of every generation.

        Parameters
        ----------
        result_dir : str
            The directory were the results are stored.
        reference_point : list
            The reference point in the original objective space. The point
            should be worse than all Pareto points for each objective.
        objective_weights : list, optional
            The optimization weights of the objectives. Use negative values
            for minimized objectives and positive values for maximized
            objectives. When no weights are provided, all objectives are
            assumed to be minimized.

        Returns
        -------
        generations : ndarray
            The generation numbers.
        hypervolume : ndarray
            The hypervolume of the Pareto front of each generation.

        """

        self.fitness_file = os.path.join(self.result_path,
                                         result_dir,
                                         'fitness.csv',
                                         )

        self.population_file = os.path.join(self.result_path,
                                            result_dir,
                                            'population.csv',
                                            )

        self.determine_pop_gen()

        n_obj = self._fitness_generations[0].shape[1]
        reference_point = np.asarray(reference_point, dtype=float)
        if reference_point.shape != (n_obj,):
            raise ValueError(
                "The reference point should contain %i values." % n_obj)

        if objective_weights is None:
            objective_weights = -np.ones(n_obj)
        objective_weights = np.asarray(objective_weights, dtype=float)
        if objective_weights.shape != (n_obj,):
            raise ValueError(
                "The objective weights should contain %i values." % n_obj)
        if np.any(objective_weights == 0.):
            raise ValueError("Objective weights should not be zero.")

        direction = -np.sign(objective_weights)
        reference_min = reference_point * direction
        hypervolume = np.zeros(self.n_gen)

        for index, generation in enumerate(self._fitness_generations):
            generation_min = generation * direction
            pareto_front = self._nondominated(generation_min)
            hypervolume[index] = self._hypervolume_minimization(
                pareto_front, reference_min)

        generations = np.arange(1, self.n_gen + 1)

        return generations, hypervolume


class PostProcessUQ:
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

    def __init__(self, case, pol_order, samples=1, method='full'):
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
        self.samples = samples
        self.method = method.lower()
        if self.method not in ('full', 'sparse'):
            raise ValueError("Unknown PCE method '%s'." % method)

    def _pce_prefix(self, objective):
        """
        Return the result file prefix for the selected PCE method.
        """
        if self.method == 'full':
            return 'full_pce_order_%i_%s' % (self.pol_order, objective)
        return 'sparse_pce_order_%i_%s_n_samples_%i' % (
            self.pol_order, objective, self.samples)

    def _pce_file(self, result_dir, objective, suffix):
        """
        Build a path to a PCE result file.
        """
        return os.path.join(
            self.result_path,
            result_dir,
            self._pce_prefix(objective) + suffix)

    @staticmethod
    def _read_pce_summary(summary_file):
        """
        Read scalar metrics from a PCE summary text file.
        """
        summary = {}
        with open(summary_file, 'r') as file:
            for line in file:
                parts = line.split()
                if not parts:
                    continue
                label = ' '.join(parts[:-1]).strip().lower()
                summary[label] = float(parts[-1])

        return summary

    @staticmethod
    def _summary_value(summary, key, summary_file):
        """
        Return a scalar metric from a PCE summary with a clear error message.
        """
        try:
            return summary[key]
        except KeyError as exc:
            raise ValueError("PCE summary file %s misses '%s'." %
                             (summary_file, key)) from exc

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

        data = pd.read_csv(distr_file)
        if data.shape[1] < 2:
            raise ValueError("Distribution file %s must contain two columns." %
                             distr_file)

        x_val = data.iloc[:, 0].to_numpy(dtype=float)
        y_val = data.iloc[:, 1].to_numpy(dtype=float)

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

        sobol_file = self._pce_file(
            result_dir, objective, '_Sobol_indices.csv')
        data = pd.read_csv(sobol_file)

        # Remove leading and trailing spaces from column names.
        data.columns = data.columns.str.strip()

        required_columns = ['name', 'Total-order Sobol indices']
        missing = [column for column in required_columns
                   if column not in data.columns]
        if missing:
            raise ValueError("Sobol file %s misses columns: %s." %
                             (sobol_file, ', '.join(missing)))

        names = data['name'].astype(str).tolist()
        sobol = data['Total-order Sobol indices'].astype(float).tolist()

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
                                'data_pdf_%s.csv' % objective
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
                                'data_cdf_%s.csv' % objective
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

        summary_file = self._pce_file(result_dir, objective, '.txt')
        summary = self._read_pce_summary(summary_file)
        loo = self._summary_value(summary, 'loo', summary_file)

        return loo

    def get_mean_std(self, result_dir, objective):
        """

        Reads the mean and standard deviation
        from the corresponding file in the result directory.

        Parameters
        ----------
        result_dir : str
            The result directory.
        objective : str
            The name of the quantity of interest.

        Returns
        -------
        mean : float
            The mean.
        std : float
            The standard deviation.

        """

        summary_file = self._pce_file(result_dir, objective, '.txt')
        summary = self._read_pce_summary(summary_file)
        mean = self._summary_value(summary, 'mean', summary_file)
        std = self._summary_value(summary, 'std. dev.', summary_file)

        return mean, std
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

        if not result_dirs:
            raise ValueError("At least one result directory is required.")

        # Store one dictionary per result directory.
        res_dict = []
        names = []
        for result_dir in result_dirs:
            names, sobol = self.get_sobol(result_dir, objective)
            res_dict.append(dict(zip(names, sobol)))

        # Get the highest Sobol index for each parameter over all directories.
        n_samples = len(res_dict)
        max_dict = {}
        for name in names:
            sobol_res = np.zeros(n_samples)
            for j, dic in enumerate(res_dict):
                if name not in dic:
                    raise ValueError(
                        "Sobol files do not contain the same parameters.")
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
