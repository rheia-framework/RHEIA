"""
The :py:mod:`pce` module contains classes to acquire the uncertainty
characterization of the stochastic parameters, to generate the random
experiment and to construct the Polynomial Chaos Expansion.
"""

import os
import multiprocessing as mp
from functools import partial
import itertools
import numpy as np
from scipy import stats, special
import sobol

class Data:
    """

    The class includes methods to create a file to store the samples
    and to store information on the stochastic design space, exctracted
    from the :file:`design_space` and :file:`stochastic_space` files.

    Parameters
    ----------
    inputs : dict
        input dictionary with information on
        the uncertainty quantification.
    space_obj : object
        object that contains information on the test case.

    """

    def __init__(self, inputs, space_obj):
        self.inputs = inputs
        self.path = os.path.split(
            os.path.dirname(
                os.path.abspath(__file__)))[0]
        self.space_obj = space_obj
        self.path_res = None

        self.stoch_data = {}
        self.filename_samples = ''

    def create_samples_file(self):
        """
        Creating the file that saved the input samples and model outputs.
        """

        # results directory
        self.path_res = os.path.join(self.path,
                                     'RESULTS',
                                     self.space_obj.case,
                                     'UQ',
                                     self.inputs['results dir'],
                                     )

        # if results directory does not exist, create one
        if not os.path.exists(self.path_res):
            os.makedirs(self.path_res)

        # file where evaluated samples are stored
        self.filename_samples = os.path.join(self.path_res, 'samples')

        if not os.path.isfile(self.filename_samples):
            with open(self.filename_samples, "w") as file:
                for name in self.stoch_data['names'] + \
                        self.inputs['objective names']:
                    file.write('%25s ' % name)
                file.write('\n')

    def read_stoch_parameters(self, var_values=[]):
        """
        Read in the stochastic design space
        and save the information in a dictionary

        Parameters
        ----------
        var_values : list, optional
            The design variable values in case of robust optimization.
            The default is [].

        """

        # check that no design variables are defined in the design space,
        # for which no corresponding deterministic values are provided in
        # var_values (only when performing robust optimization)
        if len(var_values) != len(self.space_obj.var_dict):
            raise ValueError(
                """When performing UQ, make sure that no design variableS
                   are present in design_space. Each variable should be
                   converted into a model parameter.""")

        # add values in var_values list to the variable names defined in the
        # var dict (only when performing robust optimization)
        tmp = {}
        for i, key in enumerate(self.space_obj.var_dict):
            tmp[key] = var_values[i]

        # dictionary with the deterministic parameters and deterministic values
        # for the design variables
        det_dict = {**self.space_obj.par_dict, **tmp}

        # create lists to store the bounds, mean and deviation for each
        # uncertain parameter
        l_b, u_b, mean, deviation = [], [], [], []
        for key in self.space_obj.upar_dict:
            # based on relation between deviation and mean, generate values
            # for the bounds
            if self.space_obj.upar_dict[key][0] == 'absolute':
                l_b.append(det_dict[key] - self.space_obj.upar_dict[key][2])
                u_b.append(det_dict[key] + self.space_obj.upar_dict[key][2])
                deviation.append(self.space_obj.upar_dict[key][2])
            elif self.space_obj.upar_dict[key][0] == 'relative':
                l_b.append(det_dict[key] *
                           (1. - self.space_obj.upar_dict[key][2]))
                u_b.append(det_dict[key] *
                           (1. + self.space_obj.upar_dict[key][2]))
                deviation.append(self.space_obj.upar_dict[key][2] *
                                 det_dict[key])
            else:
                raise ValueError(""" The relation of the deviation to the mean
                                     should be 'relative' or 'absolute'.""")

            mean.append(det_dict[key])

        # the dictionary with information on the uncertain parameters
        self.stoch_data = {'names': list(self.space_obj.upar_dict.keys()),
                           'types': [x[1] for x in list(
                                     self.space_obj.upar_dict.values())],
                           'l_b': l_b,
                           'u_b': u_b,
                           'mean': mean,
                           'deviation': deviation}


class RandomExperiment(Data):
    """
    RandomExperiment objects include information on the random samples,
    such as dimension, scaled and unscaled values. In addition, the
    samples can be evaluated and the deterministic model outputs are
    stored.

    Parameters
    ----------
    my_data : object
        Data object with information on the uncertain parameters.
    objective_position : int
        Index of the quantity of interest in the list of available
        model outputs.

    """

    def __init__(self, my_data, objective_position):
        self.my_data = my_data
        self.objective_position = objective_position

        self.dimension = len(my_data.stoch_data['types'])
        self.dists = [None] * len(my_data.stoch_data['types'])
        self.polytypes = [None] * len(my_data.stoch_data['types'])
        self.size = None
        self.x_u = None
        self.x_u_scaled = None

        self.y = None

        self.x_prev = None
        self.y_prev = None
        self.n_samples = None

    def n_terms(self):
        '''
        This method sets the number of samples to 2*(p+n)!/p!n!,
        i.e the number of terms in the full PC
        expansion of order p in n random variables.

        '''
        # number of uncertain parameters
        n_par = len(self.my_data.stoch_data['mean'])

        # polynomial order
        p_order = self.my_data.inputs['pol order']
        result = 1.
        mmin = min(n_par, p_order)
        for i in range(mmin):
            result *= p_order + n_par - i
        result_terms = int(result / np.math.factorial(mmin))

        # number of samples for training the PCE
        self.n_samples = 2 * result_terms

    def read_previous_samples(self, create_only_samples):
        """
        Read the previously evaluated samples
        and store them for future PCE construction.

        Parameters
        ----------
        create_only_samples : bool
            the boolean that indicates if samples
            should only be created.

        """

        with open(self.my_data.filename_samples, 'r') as file:
            lines = file.readlines()
            x_int = np.zeros((len(lines) - 1,
                              len(self.my_data.stoch_data['mean'])))
            y_int = np.zeros((len(lines) - 1, 1))

            # if samples are present in the samples file
            if len(lines) > 1:
                if (len(lines[-1].split()) !=
                    len(self.my_data.stoch_data['mean']) +
                        len(self.my_data.inputs['objective names'])):
                    raise SyntaxError(
                        """The samples file is not properly formatted
                           or it already contains samples
                           without model output.""")

                if create_only_samples:
                    raise ValueError(
                        """The samples file already contains samples
                           with model output.
                           Consider changing the result directory
                           or switching "create only samples" to False.""")

                # read in samples and deterministic model output
                for i, line in enumerate(lines[1:]):
                    line_split = line.split()
                    x_int[i] = [float(el) for el in line_split[:len(
                        self.my_data.stoch_data['mean'])]]
                    y_int[i] = float(
                        line_split[len(self.my_data.stoch_data['mean']) +
                                   self.objective_position])

        # store samples and deterministic output
        self.y_prev = np.array(y_int)
        self.x_prev = np.array(x_int)

    def create_distributions(self):
        """
        Create the distributions, polynomial distributions and polynomial types
        based on the stochastic design space. The available distributions
        are Uniform and Gaussian distributions.

        """

        polydists = [None] * len(self.my_data.stoch_data['types'])

        # iterate through the uncertain parameters defined in stochastic space
        for i, j in enumerate(self.my_data.stoch_data['types']):
            if j == 'Uniform':
                self.dists[i] = stats.uniform(
                    self.my_data.stoch_data['mean'][i] -
                    self.my_data.stoch_data['deviation'][i],
                    2. * self.my_data.stoch_data['deviation'][i])

                # scaled distribution for polynomial
                polydists[i] = stats.uniform(-1., 2.)

                # Legendre polynomial for Uniform distributions
                self.polytypes[i] = 'Legendre'

            elif j == 'Gaussian':
                self.dists[i] = stats.norm(
                    self.my_data.stoch_data['mean'][i],
                    self.my_data.stoch_data['deviation'][i])

                # scaled distribution for polynomial
                polydists[i] = stats.norm(0., 1.)

                # Legendre polynomial for Uniform distributions
                self.polytypes[i] = 'Hermite'

            else:
                raise ValueError(""" Distribution for %s not found.
                                     Choose between Uniform and Gaussian."""
                                 % self.my_data.stoch_data['types'][i])

    def create_samples(self, size=0):
        """
        Generate the samples for model evaluations. The sampling
        methods available are random sampling and Sobol sampling.
        The number of samples generated is equal to the integer
        `size`. If `size` is equal or lower than zero, no new samples
        are generated. Instead, the required samples are extracted
        from the existing samples file.

        Parameters
        ----------
        size : int, optional
            The number of newly created samples. The default is 0.

        """
        l_b = self.my_data.stoch_data['l_b']
        u_b = self.my_data.stoch_data['u_b']
        method = self.my_data.inputs['sampling method']

        # when samples need to be created
        if size > 0:
            self.x_u = np.zeros((size, self.dimension))
            if method == 'RANDOM':
                # random generation
                for i in range(self.dimension):
                    self.x_u[:, i] = self.dists[i].rvs(size)

            elif method == 'SOBOL':
                # sobol sequence
                skip = 123454
                x_tr = sobol.sample(dimension=self.dimension,
                                    n_points=size + len(self.x_prev),
                                    skip=skip)

                for i in range(self.dimension):
                    self.x_u[:, i] = self.dists[i].ppf(
                        x_tr[len(self.x_prev):, i])

            if len(self.x_prev) > 0:
                # when new samples are combined with previous samples
                self.x_u = np.concatenate((self.x_prev, self.x_u))

            self.size = (len(self.x_u))

        elif size == 0:
            # when no new samples are needed, get the ones from samples file
            self.x_u = self.x_prev
            self.size = len(self.x_u)

        else:
            # less samples needed than available in samples file
            self.x_u = np.split(self.x_prev,
                                [len(self.x_prev) + size,
                                 len(self.x_prev)])[0]
            self.size = len(self.x_u)

        # generate the scaled samples
        self.x_u_scaled = np.zeros((self.size, self.dimension))
        for i in range(self.dimension):
            for j in range(self.size):
                self.x_u_scaled[j, i] = ((self.x_u[j, i] - l_b[i]) /
                                         (u_b[i] - l_b[i]) * 2. - 1.)

    def create_only_samples(self, create_only_samples):
        """

        Add the generated samples to the samples file.
        These samples can be used for model evaluation,
        when the model is not connected to the UQ algorithm.

        Parameters
        ----------
        create_only_samples : bool
            the boolean that indicates
            if samples should only be created

        """

        if create_only_samples:
            # add the samples to the samples file, no output results
            with open(self.my_data.filename_samples, 'a+') as file:
                for i, x_u_sample in enumerate(self.x_u):
                    line = list(x_u_sample)
                    for j in line:
                        file.write('%25f ' % j)

                    file.write('\n')

    def evaluate(self, eval_func, params):
        """
        Evaluate the samples in the model
        and store the samples and outputs in the samples file.

        """

        # number of samples needed
        size = self.n_samples - len(self.x_prev)

        # create Uniform/gaussian distributions and corresponding orthogonal
        # polynomials
        self.create_samples(size=size)

        samples = self.x_u[-size:]
        unc_samples = np.tile(
            list(self.my_data.space_obj.par_dict.values()), (size, 1))

        # list with names from variables and parameters
        par_var_list = list(self.my_data.space_obj.par_dict.keys())
        par_var_list += list(self.my_data.space_obj.var_dict.keys())
        for j, elem in enumerate(self.my_data.space_obj.upar_dict.keys()):
            # get name from uncertain parameters out of upar_dict, and get the
            # index for each name out of the complete par_var_list

            unc_idx = par_var_list.index(elem)

            # based on the retrieved index, replace the deterministic value
            # of the stochastic parameter with the random sample
            unc_samples[:, unc_idx] = samples[:, j]

        eval_dict = []
        for sample in unc_samples:

            # bring parameter and variable names and values for each sample
            # into a dictionary
            eval_dict.append(
                self.my_data.space_obj.convert_into_dictionary(sample))

        if self.my_data.inputs['n jobs'] == 1:
            # linear processing
            res = []
            for index, sample in enumerate(eval_dict):
                res.append(eval_func((index + len(self.x_prev), sample), params=params))
                with open(self.my_data.filename_samples, 'a+') as file:
                    line = list(samples[index]) + res[-1]
                    for j in line:
                        file.write('%25f ' % j)
                    file.write('\n')


        else:
            # multiprocessing
            pool = mp.Pool(processes=self.my_data.inputs['n jobs'])
            res = pool.map(
                partial(
                    eval_func,
                    params=params),
                enumerate(eval_dict))
            pool.close()

            # append new samples and model outputs to samples file
            with open(self.my_data.filename_samples, 'a+') as file:
                for i, sample in enumerate(samples):
                    line = list(np.concatenate((sample, res[i])))
                    for j in line:
                        file.write('%25f ' % j)

                    file.write('\n')

        # check that the quantity of interest exists
        if self.objective_position > len(res[0]) - 1:
            raise IndexError(""" The objective "%s" falls out of
                                 the range of predefined quantities
                                 of interest. Only %i outputs are
                                 returned from the model""" %
                             (self.my_data.inputs['objective names'][
                                 int(self.objective_position)], len(res[0])))

        # combine model output for quantity of interest with previous results
        # available in samples file
        y_res = [row[self.objective_position] for row in res]
        if self.y_prev.size:
            self.y = np.vstack((self.y_prev, np.array(y_res).reshape((-1, 1))))
        else:
            self.y = np.array(y_res).reshape((-1, 1))

class PCE(RandomExperiment):
    """

    Class which creates a Polynomial Chaos Expansion (PCE) object. A PCE is
    characterized by the following attributes:

        - Basis : PC basis functions
        - Coefficients : PC coefficients
        - Moments : Statistical moments
        - sensitivity : Sobol' indices


    Parameters
    ----------
    RandomExperiment : obj
        RandomExperiment object, with information on the random samples.

    """

    def __init__(self, Experiment):
        self.my_experiment = Experiment
        self.order = self.my_experiment.my_data.inputs['pol order']
        self.a_matrix = 0.0
        self.loo = 0.0
        self.basis = dict()
        self.coefficients = dict()
        self.sensitivity = dict()
        self.moments = dict()
        self.psi_sq = 0.0

    #####################################
    # module uncertainty quantification #
    #####################################

    def n_to_sum(self, n, s):
        '''
        This function creates a list
        of all possible vectors of length 'n' that sum to 's'

        Parameters
        ----------
        n : int
            vector's length (dimension)
        s : int
            sum (total order)

        Returns
        -------
        The output is a generator object to be called as follows
            - list(n_to_sum(dimension,i))
        '''

        if n == 1:
            yield (s,)
        else:
            for i in range(s + 1):
                for j in self.n_to_sum(n - 1, s - i):
                    yield j + (i,)

    def multindices(self, idx):
        """

        This method returns a set of multi-indices

        Parameters
        ----------
        idx : list
            the range for the number of terms in the PCE

        Returns
        -------
        list
            the list with multi-indices

        """
        dimension = self.my_experiment.dimension

        multindices = list()

        # generate the set of multi-indices
        for i in range(self.order + 1):
            multindices.extend(list(self.n_to_sum(dimension, i)))

        return [multindices[i] for i in idx]

    def ols(self, a_matrix, b_matrix):
        """

        Perform Ordinary Least Squares on the input matrices a_matrix and
        b_matrix.
        Return the result.

        Parameters
        ----------
        a_matrix : array
            the matrix with the response of the basis functions to the samples
        b_matrix : array
            the model output for the input samples

        Returns
        -------
        2_D array
            result of the ordinary least-square regression

        """

        m, n = a_matrix.shape
        deg_of_freedom = m - n

        sol = np.linalg.lstsq(a_matrix, b_matrix, rcond=None)
        sol_2 = np.column_stack(1.0 / deg_of_freedom * sol[1])

        var = np.dot(np.row_stack(np.diag(np.linalg.inv(
            np.dot(np.transpose(a_matrix), a_matrix)))), sol_2)

        return np.column_stack((sol[0], var))

    def calc_a(self, multindices):
        """

        This method builds the matrix containing the basis functions evaluated
        at sample locations, i.e. the matrix A in Au = b

        Parameters
        ----------
        multindices : list
            the list with the multi-indices

        Returns
        -------
        A : array
            the matrix with the evaluated basis functions for the samples

        """

        dimension = self.my_experiment.dimension
        x_u_scaled = self.my_experiment.x_u_scaled
        a_matrix = np.ones([self.my_experiment.size, len(multindices)])

        # generate the A matrix
        for i, multiindex in enumerate(multindices):
            for j in range(dimension):
                deg = multiindex[j]
                if self.my_experiment.polytypes[j] == 'Legendre':
                    a_matrix[:, i] *= special.eval_legendre(deg,
                                                            x_u_scaled[:, j])

                elif self.my_experiment.polytypes[j] == 'Hermite':
                    a_matrix[:, i] *= special.eval_hermitenorm(deg,
                                                               x_u_scaled[:,
                                                                          j])

        return a_matrix

    def run(self):
        """

        Solve Ordinary Least Squares problem
        Full PC expansion is assumed containing n_terms(dimension,order)


        """

        # number of terms in the PCE
        n_terms = int(self.my_experiment.n_samples / 2)

        self.basis['model'] = range(n_terms)

        # full PCE basis
        self.basis['multi-indices'] = self.multindices(self.basis['model'])
        self.basis['polytypes'] = self.my_experiment.polytypes

        # determine the A matrix
        self.a_matrix = self.calc_a(self.basis['multi-indices'])

        # get the coefficient through OLS
        self.coefficients = self.ols(self.a_matrix, self.my_experiment.y)

        # acquire the statistics
        self.get_statistics(mean=True, variance=True)

    ################################
    # module statistics extraction #
    ################################

    def get_statistics(self, mean=True, variance=True):
        '''
        This function calculates high order moments (up to order 2) by
        taking advantage of the fact that any permutation of indices
        will lead to the same value for the summand.

        Returns Statistics: mean, std. deviation

        '''

        coeff = self.coefficients

        # quantify the mean
        if mean:
            if self.basis['model'][0] == 0:
                self.moments['mean'] = coeff[0, 0]

        if variance:
            # the term <psi_i, psi_j>
            self.psi_sq = self.get_psi_sq()

            # quantify the variance
            var = 0.0
            if self.basis['model'][0] == 0:
                for i in range(1, len(self.basis['model'])):
                    var += self.psi_sq[i] * coeff[i, 0]**2

            self.moments['variance'] = var

    def get_psi_sq(self):
        """

        Calculate the term <psii,psij>

        Returns
        -------
        psi_sq : array
            the term <psii,psij>

        """
        dim = self.my_experiment.dimension

        multindices = self.basis['multi-indices']

        n_terms = len(multindices)

        psi_sq = np.ones(n_terms,)

        for i in range(n_terms):
            for j in range(dim):
                deg = multindices[i][j]

                if self.my_experiment.polytypes[j] == 'Legendre':

                    x_i, w_i = special.p_roots(deg + 1)

                    '''
                    Integrate exactly the SQUARE
                    of the Legendre polynomial. For example,
                    if the Legendre polynomial is of order (deg),
                    the numerical integration must be exact
                    till order (deg**2). Thus, we need at least
                    (deg+1) abscissas' and weights.
                    '''
                    poly = special.legendre(deg)**2
                    psi_sq[i] *= 1.0 / 2 * sum(w_i * poly(x_i))

                elif self.my_experiment.polytypes[j] == 'Hermite':

                    x_i, w_i = special.he_roots(deg + 1)

                    '''
                    special.he_roots(deg) and
                    np.polynomial.hermite_e.hermegauss(deg)
                    returns the same abscissas'
                    but different weights (!). There is a factor 2
                    between the two. Given the fact that the integral of
                    the standard Gaussian must be 1,
                    np.polynomial.hermite_e.hermegauss(deg)
                    provides the right weights.
                    '''
                    poly = special.hermitenorm(deg)**2

                    # 2*w_i*poly(x_i)
                    psi_sq[i] *= 1.0 / np.sqrt(2 * np.pi) * sum(w_i *
                                                                poly(x_i))

        return psi_sq

    ######################
    # module sensitivity #
    ######################

    def calc_sobol(self):
        '''
        This method calculates the Sobol' indices Si of PCE.
        They represent the fraction of the total variance
        that can be attributed to each input variable i (Si)
        or combinations of variables (i1,i2,...,is).
        The total Sobol indices quantify the total effect
        of variable i, i.e. including the first order effect
        and the interactions with the other variables.

        '''

        dimension = self.my_experiment.dimension

        var = self.moments['variance']

        coeff = self.coefficients

        psi_sq = self.get_psi_sq()  # vector of size (P+1)

        # Create list of all of terms
        terms = list()
        for i in range(1, self.order + 1):
            terms.extend(list(itertools.combinations(range(dimension), i)))

        copy = self.basis['multi-indices'][:]
        del copy[0]

        # Compute the Sobol' indices
        s_i = np.zeros(len(terms),)
        for i, multindex in enumerate(copy):
            for j, term in enumerate(terms):
                if (len(np.nonzero(multindex)[0]) == len(term)) and (
                        np.nonzero(multindex)[0] == term).all():
                    s_i[j] += 1.0 / var * coeff[i + 1][0]**2 * psi_sq[i + 1]

                    break

        # Compute the total Sobol' indices
        s_tot_i = np.zeros(dimension,)
        for i in range(dimension):
            idx = np.array([item for item in range(len(terms))
                            if (terms[item] == np.array([i])).any()])
            for index, j in enumerate(idx):
                s_tot_i[i] += s_i[j]

        self.sensitivity['s_i'] = s_i
        self.sensitivity['s_tot_i'] = s_tot_i
        self.sensitivity['terms'] = terms

    ##############################
    # accuracy evaluation module #
    ##############################

    def calc_loo(self):
        '''
        This method evaluates the Leave-One-Out (LOO) error for the
        constructed PCE. The formula is adopted from Sudret et al. [1].

        [1] Sudret, B. (2014). Polynomial chaos expansions and stochastic
        finite-element methods. In Risk and Reliability in Geotechnical
        Engineering.

        '''

        size = self.my_experiment.size
        y_res = self.my_experiment.y
        y_hat = np.row_stack((np.dot(self.a_matrix, self.coefficients[:, 0])))

        # acquire the diagonal matrix D
        b_matrix = np.linalg.inv(np.dot(np.transpose(self.a_matrix),
                                        self.a_matrix))
        c_matrix = np.dot(b_matrix, np.transpose(self.a_matrix))
        h_matrix = np.dot(self.a_matrix, c_matrix)
        d_matrix = np.diag(h_matrix)

        # quantify the LOO error
        self.loo = 0.
        for i in range(size):
            deltai = (y_res[i, 0] - y_hat[i, 0]) / (1 - d_matrix[i])

            self.loo += 1.0 / float(size) * deltai**2. / (np.var(y_res))

        if self.loo > 1.:
            raise Warning("""The LOO error is higher than 1.
                                Check the UQ characterization and results.""")

    ##########################
    # result printing module #
    ##########################

    def print_res(self):
        """
        This method prints an overview of the inputs and results of the PCE.
        Additionally, it writes the PCE information and Sobol' indices
        in the result files.'
        """

        mean = self.moments['mean']
        var = self.moments['variance']

        print('%s Polynomial chaos output %s \n' % ('-' * 20, '-' * 20))
        print(' Number of input variables ::'.ljust(30) + '%d' %
              len(self.my_experiment.dists))
        print(' Maximal degree ::'.ljust(30) + '%d' %
              max([sum(i) for i in self.basis['multi-indices']]))
        print(' Size of full basis ::'.ljust(30) +
              '%d' % self.my_experiment.n_samples)
        print(' Size of sparse basis ::'.ljust(30) +
              '%d' % len(self.basis['multi-indices']) + '\n')

        print(' Full model evaluation ::'.ljust(30) +
              '%d' % self.my_experiment.size)
        print(' Leave-one-out error ::'.ljust(30) + '%s' % self.loo + '\n')

        print(' Mean value ::'.ljust(30) + '%s' % mean)
        print(' Std. deviation ::'.ljust(30) + '%s' % np.sqrt(var))
        print(' First-order Sobol indices ::'.ljust(30) +
              '%s' % np.around(self.sensitivity['s_i'], decimals=4) + '\n')
        print(' Total Sobol indices ::'.ljust(30) +
              '%s' % np.around(self.sensitivity['s_tot_i'], decimals=4) + '\n')

        print('-' * 65)

        filename_res = (
            "full_pce_order_%d_%s" %
            (self.order,
             self.my_experiment.my_data.inputs['objective of interest']))

        # write statistics in the result file
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               filename_res), "w") as file:
            file.write('%25s %25f \n' % ('LOO', self.loo))
            file.write(
                '%25s %25f \n' %
                ('sparse basis', len(
                    self.basis['multi-indices'])))
            file.write('%25s %25f \n' % ('mean', mean))
            file.write(
                '%25s %25f \n' % ('std. dev.', np.sqrt(var)))

        # write sobol indices in the corresponding result file
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               filename_res + '_Sobol_indices'), "w") as file:
            file.write(
                '%30s %30s %30s  \n' %
                ('name',
                 'First-order Sobol indices',
                 'Total-order Sobol indices'))
            indices = np.argsort(self.sensitivity['s_tot_i'])[::-1]
            for i in indices:
                file.write(
                    '%30s %30f %30f  \n' %
                    (self.my_experiment.my_data.stoch_data['names'][i],
                     self.sensitivity['s_i'][i],
                     self.sensitivity['s_tot_i'][i]))

    def draw(self, size):
        """
        This module creates the probability density function
        and cumulative distribution function and writes
        the corresponding values to construct these values
        in the result files.

        Parameters
        ----------
        size : int
            number of samples to be created

        """

        # create samples to evaluate on the PCE
        self.my_experiment.create_samples(size=size)

        # generate A matrix and the coefficients u
        a_matrix = self.calc_a(self.basis['multi-indices'])
        x_u_scaled = np.row_stack((self.coefficients[:, 0]))

        data = np.dot(a_matrix, x_u_scaled)

        # generate the pdf
        density, bins = np.histogram(data, bins=100, density=1)
        centers = [(a + b) / 2 for a, b in zip(bins[::1], bins[1::1])]
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               ("data_pdf_%s"
                                % (self.my_experiment.my_data.inputs[
                                    'objective of interest']))), "w") as file:

            file.write(
                '%25s %25s \n' %
                (self.my_experiment.my_data.inputs['objective of interest'],
                 'probability density'))
            for i, j in enumerate(centers):
                file.write('%25f %25f' % (j, density[i]))
                file.write('\n')

        # generate the cdf
        cdf = np.cumsum(density * np.diff(bins))
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               ("data_cdf_%s"
                                % (self.my_experiment.my_data.inputs[
                                    'objective of interest']))), "w") as file:

            file.write(
                '%25s %25s \n' %
                (self.my_experiment.my_data.inputs['objective of interest'],
                 'cumulative probability'))
            for i, j in enumerate(centers):
                file.write('%25f %25f' % (j, cdf[i]))
                file.write('\n')
