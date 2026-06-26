"""
The :py:mod:`pce` module contains classes to acquire the uncertainty
characterization of the stochastic parameters, to generate the random
experiment and to construct the Polynomial Chaos Expansion.
"""

import os
import csv
import multiprocessing as mp
from functools import partial
import itertools
import numpy as np
from scipy import stats, special
from scipy.stats import qmc
import warnings
import pandas as pd
import math


def _as_float_scalar(value, name):
    """
    Convert a scalar-like model value to float.
    """
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError("%s should contain exactly one scalar value." % name)

    return float(array.reshape(-1)[0].item())


def _normalize_outputs(results):
    """
    Normalize model outputs to a rectangular list of floats.
    """
    normalized = []
    for row_index, row in enumerate(results):
        if np.isscalar(row) or np.asarray(row).ndim == 0:
            normalized.append([
                _as_float_scalar(row, "Model output %i" % row_index)])
            continue

        try:
            values = list(row)
        except TypeError as exc:
            raise TypeError(
                "Model output %i should be a scalar or a sequence." %
                row_index) from exc

        normalized.append([
            _as_float_scalar(value, "Model output %i item %i" %
                             (row_index, value_index))
            for value_index, value in enumerate(values)])

    if not normalized:
        raise ValueError("No model outputs were returned.")

    n_outputs = len(normalized[0])
    if n_outputs == 0:
        raise ValueError("Model outputs should not be empty.")
    if any(len(row) != n_outputs for row in normalized):
        raise ValueError("All model output rows should have the same length.")

    return normalized


def _lognormal_to_normal(mean, deviation):
    """
    Return latent normal parameters from lognormal mean and standard deviation.
    """
    if mean <= 0.:
        raise ValueError("A lognormal distribution requires a positive mean.")
    if deviation <= 0.:
        raise ValueError(
            "A lognormal distribution requires a positive standard deviation.")

    sigma = np.sqrt(np.log(1. + (deviation / mean) ** 2))
    mu = np.log(mean) - 0.5 * sigma ** 2

    return mu, sigma


class Data:
    """

    The class includes methods to create a file to store the samples
    and to store information on the stochastic design space, exctracted
    from the :file:`design_space.txt` and :file:`stochastic_space.txt` files.

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
        Create the result directory and samples file for the UQ run.
        """

        self.path_res = os.path.join(self.path,
                                     'RESULTS',
                                     self.space_obj.case,
                                     'UQ',
                                     self.inputs['results dir'],
                                     )

        os.makedirs(self.path_res, exist_ok=True)

        self.filename_samples = os.path.join(self.path_res, 'samples.csv')
        
        if not os.path.isfile(self.filename_samples):
            columns = self.stoch_data['names'] + self.inputs['objective names']
            df = pd.DataFrame([columns])
            with open(self.filename_samples, 'w') as f:
                df.to_csv(f, header=False, index=False, lineterminator='\n')
        
        
    def read_stoch_parameters(self, var_values=None):
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
        if var_values is None:
            var_values = []

        if len(var_values) != len(self.space_obj.var_dict):
            raise ValueError(
                """When performing UQ, make sure that no design variableS
                   are present in design_space. Each variable should be
                   converted into a model parameter.""")

        # add values in var_values list to the variable names defined in the
        # var dict (only when performing robust optimization)
        tmp = {key: var_values[i] for i, key in enumerate(self.space_obj.var_dict)}

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
        self._lognorm_mu = [None] * len(my_data.stoch_data['types'])
        self._lognorm_sigma = [None] * len(my_data.stoch_data['types'])
        self.size = None
        self.x_u = None
        self.x_u_scaled = None

        self.y = None

        self.x_prev = None
        self.y_prev = None
        self.n_samples = None
        self.n_terms_full = None

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
        result_terms = int(result / math.factorial(mmin))
        self.n_terms_full = result_terms

        uq_method = self.my_data.inputs.get('uq method', 'full')
        self.my_data.inputs['uq method'] = uq_method

        if uq_method == 'full':
            # number of samples for training the PCE
            self.n_samples = 2 * result_terms
        elif uq_method == 'sparse':
            self.n_samples = int(self.my_data.inputs['n samples'])
        else:
            raise ValueError("""The UQ method should be 'full' or 'sparse'.""")


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
        
        n_inputs = len(self.my_data.stoch_data['mean'])
        n_outputs = len(self.my_data.inputs['objective names'])

        with open(self.my_data.filename_samples, 'r') as file:
            lines = file.readlines()

        previous_lines = [line.strip().split(",") for line in lines[1:]]
        n_previous = len(previous_lines)

        self.x_prev = np.zeros((n_previous, n_inputs))
        self.y_prev = np.zeros((n_previous, 1))

        if n_previous == 0:
            return

        row_lengths = {len(line) for line in previous_lines}
        input_only = row_lengths == {n_inputs}
        evaluated = row_lengths == {n_inputs + n_outputs}

        if not input_only and not evaluated:
            raise SyntaxError(
                """The samples file is not properly formatted.""")

        if create_only_samples and evaluated:
            raise ValueError(
                """The samples file already contains samples with model output.
                   Consider changing the result directory or switching
                   "create only samples" to False.""")

        if not create_only_samples and input_only:
            raise SyntaxError(
                """The samples file already contains samples without model output.
                   Evaluate these samples first or use a different result
                   directory.""")

        for i, line_split in enumerate(previous_lines):
            self.x_prev[i] = [float(el) for el in line_split[:n_inputs]]
            if evaluated:
                self.y_prev[i] = float(
                    line_split[n_inputs + self.objective_position])
            
    def create_distributions(self):
        """
        Create the distributions, polynomial distributions and polynomial types
        based on the stochastic design space. The available distributions
        are Uniform, Gaussian and Lognormal distributions.

        """

        # iterate through the uncertain parameters defined in stochastic space
        for i, j in enumerate(self.my_data.stoch_data['types']):
            if j == 'Uniform':
                self.dists[i] = stats.uniform(
                    self.my_data.stoch_data['mean'][i] -
                    self.my_data.stoch_data['deviation'][i],
                    2. * self.my_data.stoch_data['deviation'][i])

                # Legendre polynomial for Uniform distributions
                self.polytypes[i] = 'Legendre'

            elif j == 'Gaussian':
                self.dists[i] = stats.norm(
                    self.my_data.stoch_data['mean'][i],
                    self.my_data.stoch_data['deviation'][i])

                # Hermite polynomial for Gaussian distributions
                self.polytypes[i] = 'Hermite'

            elif j == 'Lognormal':
                mu, sigma = _lognormal_to_normal(
                    self.my_data.stoch_data['mean'][i],
                    self.my_data.stoch_data['deviation'][i])
                self._lognorm_mu[i] = mu
                self._lognorm_sigma[i] = sigma
                self.dists[i] = stats.lognorm(s=sigma, scale=np.exp(mu))

                # Hermite polynomial for the latent Gaussian distribution
                self.polytypes[i] = 'Hermite'

            else:
                raise ValueError(""" Distribution for %s not found.
                                     Choose between Uniform, Gaussian and
                                     Lognormal."""
                                 % self.my_data.stoch_data['types'][i])

    def _scale_samples_for_pce(self):
        """
        Transform physical samples to the domain used by the PCE basis.
        """
        l_b = np.asarray(self.my_data.stoch_data['l_b'], dtype=float)
        u_b = np.asarray(self.my_data.stoch_data['u_b'], dtype=float)
        x_u_scaled = np.zeros_like(self.x_u, dtype=float)

        for i, dist_type in enumerate(self.my_data.stoch_data['types']):
            if dist_type == 'Lognormal':
                if np.any(self.x_u[:, i] <= 0.):
                    raise ValueError(
                        "Lognormal samples should be strictly positive.")
                mu = self._lognorm_mu[i]
                sigma = self._lognorm_sigma[i]
                if mu is None or sigma is None:
                    mu, sigma = _lognormal_to_normal(
                        self.my_data.stoch_data['mean'][i],
                        self.my_data.stoch_data['deviation'][i])
                    self._lognorm_mu[i] = mu
                    self._lognorm_sigma[i] = sigma
                x_u_scaled[:, i] = (np.log(self.x_u[:, i]) - mu) / sigma
            else:
                denom = u_b[i] - l_b[i]
                if denom == 0.:
                    raise ValueError(
                        "At least one stochastic parameter has u_b == l_b.")
                x_u_scaled[:, i] = (self.x_u[:, i] - l_b[i]) / denom * 2. - 1.

        return x_u_scaled

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
        method = self.my_data.inputs['sampling method']

        size = int(size)

        if size > 0:
            self.x_u = np.zeros((size, self.dimension))
            if method == 'RANDOM':
                for i in range(self.dimension):
                    self.x_u[:, i] = self.dists[i].rvs(size)

            elif method == 'SOBOL':
                sampler = qmc.Sobol(
                    d=self.dimension,
                    scramble=True,
                    rng=42,
                )

                n_previous = len(self.x_prev)

                if n_previous > 0:
                    sampler.fast_forward(n_previous)
        
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message="The balance properties of Sobol' points")
                    x_tr = sampler.random(n=size)

                for i in range(self.dimension):
                    self.x_u[:, i] = self.dists[i].ppf(x_tr[:, i])

            else:
                raise ValueError("""The sampling method should be 'RANDOM'
                                     or 'SOBOL'.""")

            if len(self.x_prev) > 0:
                self.x_u = np.concatenate((self.x_prev, self.x_u))

            self.size = len(self.x_u)

        elif size == 0:
            self.x_u = self.x_prev
            self.size = len(self.x_u)

        else:
            self.x_u = self.x_prev[:len(self.x_prev) + size]
            self.size = len(self.x_u)

        # Map physical samples to the domain used by the polynomial basis.
        self.x_u_scaled = self._scale_samples_for_pce()

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

        if create_only_samples and self.size > 0:
            n_previous = len(self.x_prev)
            samples = self.x_u[n_previous:]
            if len(samples) == 0:
                return

            df = pd.DataFrame(samples, columns=None)
            with open(self.my_data.filename_samples, 'a+') as f:
                df.to_csv(f, header=False, index=False, lineterminator='\n')

    def evaluate(self, eval_func, params):
        """
        Evaluate the samples in the model
        and store the samples and outputs in the samples file.

        """

        size = self.n_samples - len(self.x_prev)
        if size <= 0:
            return

        # Create only the new input samples. Previous samples remain in x_prev.
        self.create_samples(size=size)

        samples = self.x_u[-size:]
        unc_samples = np.tile(
            list(self.my_data.space_obj.par_dict.values()), (size, 1))

        par_var_list = list(self.my_data.space_obj.par_dict.keys())
        par_var_list += list(self.my_data.space_obj.var_dict.keys())

        for j, elem in enumerate(self.my_data.space_obj.upar_dict.keys()):
            unc_idx = par_var_list.index(elem)
            unc_samples[:, unc_idx] = samples[:, j]

        eval_dict = [
            self.my_data.space_obj.convert_into_dictionary(sample)
            for sample in unc_samples]

        expected_n_outputs = len(self.my_data.inputs['objective names'])

        def normalize_result(result):
            normalized = _normalize_outputs([result])[0]

            if len(normalized) != expected_n_outputs:
                raise ValueError(
                    "The model returned %i outputs, but %i objective names "
                    "were provided in the run dictionary." %
                    (len(normalized), expected_n_outputs))

            if self.objective_position > len(normalized) - 1:
                raise IndexError(""" The objective "%s" falls out of
                                     the range of predefined quantities
                                     of interest. Only %i outputs are
                                     returned from the model""" %
                                 (self.my_data.inputs['objective names'][
                                     int(self.objective_position)],
                                  len(normalized)))

            return normalized

        y_res = []
        with open(self.my_data.filename_samples, 'a+', newline='') as file:
            writer = csv.writer(file, lineterminator='\n')

            if self.my_data.inputs['n jobs'] == 1:
                result_iter = (
                    eval_func((index + len(self.x_prev), sample),
                              params=params)
                    for index, sample in enumerate(eval_dict))

                for sample, result in zip(samples, result_iter):
                    normalized = normalize_result(result)
                    writer.writerow(list(sample) + normalized)
                    file.flush()
                    y_res.append(normalized[self.objective_position])

            else:
                with mp.Pool(processes=self.my_data.inputs['n jobs']) as pool:
                    result_iter = pool.imap(
                        partial(eval_func, params=params),
                        enumerate(eval_dict))

                    for sample, result in zip(samples, result_iter):
                        normalized = normalize_result(result)
                        writer.writerow(list(sample) + normalized)
                        file.flush()
                        y_res.append(normalized[self.objective_position])

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

    def forward(self, model, candidates, multindices):
        # t = time.time()
        T = np.zeros(len(candidates), )

        '''
        #if parallel
        results = self.uq_multiprocessing(candidates,multindices)      
        T = results[:]
        '''
        # if not parallel
        n = self.my_experiment.size
        dof = n - 1

        for i, candidate in enumerate(candidates):
            A = np.zeros([n, 1])
            A[:, 0] = self.a_matrix[:, candidate]

            dummy = np.linalg.inv(np.dot(np.transpose(A), A))
            beta_OLS = np.dot(np.dot(dummy, np.transpose(A)), self.residu)

            rhat = beta_OLS * A

            s2 = 1.0 / dof * np.sum((self.residu[:, 0] - rhat[:, 0]) ** 2)
            var_OLS = dummy * s2

            T[i] = beta_OLS.item() / np.sqrt(var_OLS.item())
            ###
        # '''
        jmax = np.argmax(np.abs(T))

        #        print '     --- Term n°%d was added to the model %s'%(candidates[jmax], multindices[candidates[jmax]])

        model.append(candidates[jmax])

        candidates = list(candidates)  # Convert range to a list
        candidates.remove(candidates[jmax])
        # print time.time() - t

        return model, candidates

    def backward(self, model, candidates, multindices):
        alpha = 0.05

        n = self.my_experiment.size
        p = len(model)

        tstar = stats.t.ppf(1 - alpha / 2, n - p)

        C = self.coefficients

        A = np.sign(C[:, 0] - tstar * np.sqrt(C[:, 1])) * np.sign(C[:, 0] + tstar * np.sqrt(C[:, 1]))

        ind2remove = np.where(A < 0)[0]
        if len(ind2remove) != 0:
            #       for i in range(len(ind2remove)):
            #                print '     --- Term n°%d was removed from the model %s'%(model[ind2remove[i]], multindices[model[ind2remove[i]]])

            ind2remove = [
                index for index in ind2remove
                if index < len(model) and model[index] != 0]
            model = [x for i, x in enumerate(model) if i not in ind2remove]

        return model, candidates

    def truncate(self, model, coeff, nb):

        if 0 in model:
            intercept_pos = model.index(0)
            keep_count = max(nb - 1, 0)
        else:
            intercept_pos = None
            keep_count = nb

        A = np.divide(
            np.abs(np.sqrt(coeff[:, 1])),
            np.abs(coeff[:, 0]),
            out=np.full(len(coeff), np.inf),
            where=coeff[:, 0] != 0.) * 100
        if intercept_pos is not None:
            A[intercept_pos] = np.inf

        sortedA = np.argsort(A)

        output = []
        if intercept_pos is not None:
            output.append(0)
        output.extend(model[i] for i in sortedA[:keep_count])

        return output


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
        if deg_of_freedom <= 0:
            raise ValueError(
                "OLS needs more samples than basis terms. "
                "Received %i samples and %i basis terms." % (m, n))

        sol = np.linalg.lstsq(a_matrix, b_matrix, rcond=None)
        residual = sol[1]
        if residual.size == 0:
            residual = np.sum((b_matrix - np.dot(a_matrix, sol[0])) ** 2,
                              axis=0)
        sol_2 = np.atleast_2d(1.0 / deg_of_freedom * residual)

        var = np.dot(np.vstack(np.diag(np.linalg.pinv(
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

    def run(self, uq_method=None):
        """
        Fit a full or sparse PCE to the evaluated random experiment.

        """

        # number of terms in the full PCE basis
        if uq_method is None:
            uq_method = self.my_experiment.my_data.inputs.get('uq method', 'full')
        self.my_experiment.my_data.inputs['uq method'] = uq_method

        if self.my_experiment.n_terms_full is None:
            self.my_experiment.n_terms()

        n_terms = int(self.my_experiment.n_terms_full)
        self.full_basis = n_terms

        self.residu = self.my_experiment.y

        if uq_method == 'full':

            self.basis['model'] = range(n_terms)

            # full PCE basis
            self.basis['multi-indices'] = self.multindices(self.basis['model'])
            self.basis['polytypes'] = self.my_experiment.polytypes

            # determine the A matrix
            self.a_matrix = self.calc_a(self.basis['multi-indices'])

            # get the coefficient through OLS
            self.coefficients = self.ols(self.a_matrix, self.my_experiment.y)

        elif uq_method == 'sparse':

            candidates = list(range(1, n_terms))  # keep intercept in model
            multindices = self.multindices(candidates)
            all_multindices = self.multindices(range(n_terms))
            full_a_matrix = self.calc_a(all_multindices)
            self.a_matrix = full_a_matrix

            niter = int(self.my_experiment.n_samples / 2)
            model = [0]
            self.coefficients = self.ols(full_a_matrix[:, model],
                                         self.my_experiment.y)
            y_hat = np.vstack(np.dot(full_a_matrix[:, model],
                                     self.coefficients[:, 0]))
            self.residu = self.my_experiment.y - y_hat

            # Forward/backward selection builds a candidate sparse model.
            for _ in range(niter):
                if not candidates:
                    break
                model, candidates = self.forward(model, candidates, multindices)
                self.coefficients = self.ols(full_a_matrix[:, model],
                                             self.my_experiment.y)
                model, candidates = self.backward(model, candidates, multindices)
                self.coefficients = self.ols(full_a_matrix[:, model],
                                             self.my_experiment.y)
                y_hat = np.vstack(np.dot(full_a_matrix[:, model],
                                         self.coefficients[:, 0]))
                self.residu = self.my_experiment.y - y_hat

            selected_model = list(model)
            selected_coefficients = self.coefficients

            best = 0
            current = 10 ** 6
            for i in range(len(selected_model)):

                model = self.truncate(selected_model, selected_coefficients, i + 1)
                model.sort()
                self.a_matrix = full_a_matrix[:, model]
                self.coefficients = self.ols(self.a_matrix, self.my_experiment.y)

                self.basis['model'] = model
                self.basis['multi-indices'] = [all_multindices[idx]
                                               for idx in model]
                self.basis['polytypes'] = self.my_experiment.polytypes

                self.calc_loo(warn=False)
                # print self.LOO
                if self.loo < current:
                    current = self.loo
                    best = i + 1
                    # break

                self.a_matrix = full_a_matrix

            model = self.truncate(selected_model, selected_coefficients, best)
            model.sort()
            self.a_matrix = full_a_matrix[:, model]

            self.coefficients = self.ols(self.a_matrix, self.my_experiment.y)
            self.basis['model'] = model
            self.basis['multi-indices'] = [all_multindices[idx]
                                           for idx in model]
            self.basis['polytypes'] = self.my_experiment.polytypes

        else:
            raise ValueError("""The UQ method should be 'full' or 'sparse'.""")

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
        constant_position = (
            self.basis['model'].index(0) if 0 in self.basis['model'] else None)

        if mean:
            self.moments['mean'] = (
                coeff[constant_position, 0]
                if constant_position is not None else 0.0)

        if variance:
            # the term <psi_i, psi_j>
            self.psi_sq = self.get_psi_sq()

            # quantify the variance
            var = 0.0
            for i in range(len(self.basis['model'])):
                if i == constant_position:
                    continue
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

        indexed_terms = [
            (index, multindex)
            for index, multindex in enumerate(self.basis['multi-indices'])
            if any(multindex)]

        if var == 0.:
            raise ValueError(
                "Sobol' indices cannot be calculated for zero variance.")

        # Compute the Sobol' indices
        s_i = np.zeros(len(terms),)
        for coeff_index, multindex in indexed_terms:
            for j, term in enumerate(terms):
                if (len(np.nonzero(multindex)[0]) == len(term)) and (
                        np.nonzero(multindex)[0] == term).all():
                    s_i[j] += (1.0 / var * coeff[coeff_index][0]**2 *
                               psi_sq[coeff_index])

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

    def calc_loo(self, warn=True):
        '''
        This method evaluates the Leave-One-Out (LOO) error for the
        constructed PCE. The formula is adopted from Sudret et al. [1].

        [1] Sudret, B. (2014). Polynomial chaos expansions and stochastic
        finite-element methods. In Risk and Reliability in Geotechnical
        Engineering.

        '''

        y_res = self.my_experiment.y
        y_hat = np.vstack((np.dot(self.a_matrix, self.coefficients[:, 0])))
        y_var = np.var(y_res)
        if y_var == 0.:
            raise ValueError("LOO error cannot be calculated for zero variance.")

        # acquire the diagonal matrix D
        b_matrix = np.linalg.pinv(np.dot(np.transpose(self.a_matrix),
                                         self.a_matrix))
        c_matrix = np.dot(b_matrix, np.transpose(self.a_matrix))
        h_matrix = np.dot(self.a_matrix, c_matrix)
        d_matrix = np.diag(h_matrix)

        # quantify the LOO error
        delta = (y_res[:, 0] - y_hat[:, 0]) / (1. - d_matrix)
        self.loo = np.mean(delta**2.) / y_var

        if warn and self.loo > 1.:
            warnings.warn("The LOO error is higher than 1. "
                          "Check the UQ characterization and results.")

    ##########################
    # result printing module #
    ##########################

    def result_filename(self):
        """
        Return the base result filename for the fitted PCE.
        """
        method = self.my_experiment.my_data.inputs['uq method']
        objective = self.my_experiment.my_data.inputs['objective of interest']

        if method == 'sparse':
            n_samples = self.my_experiment.my_data.inputs['n samples']
            return "%s_pce_order_%d_%s_n_samples_%i.txt" % (
                method, self.order, objective, n_samples)

        return "%s_pce_order_%d_%s.txt" % (
            method, self.order, objective)

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
              '%d' % self.full_basis)
        print(' Size of sparse basis ::'.ljust(30) +
              '%d' % len(self.basis['multi-indices']) + '\n')

        print(' Full model evaluation ::'.ljust(30) +
              '%d' % self.my_experiment.size)
        print(' Leave-one-out error ::'.ljust(30) + '%s' % self.loo + '\n')

        print(' Mean value ::'.ljust(30) + '%s' % mean)
        print(' Std. deviation ::'.ljust(30) + '%s' % np.sqrt(var))
        print(' First-order Sobol indices ::'.ljust(30) +
              '%s' % np.around(self.sensitivity['s_i'][
                  :len(self.sensitivity['s_tot_i'])], decimals=4) + '\n')
        print(' Total Sobol indices ::'.ljust(30) +
              '%s' % np.around(self.sensitivity['s_tot_i'], decimals=4) + '\n')

        print('-' * 65)

        filename_res = self.result_filename()

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
                               filename_res[:-4] + '_Sobol_indices.csv'), "w") as file:                               

            file.write(
                '%s,%s,%s  \n' %
                ('name',
                 'First-order Sobol indices',
                 'Total-order Sobol indices'))
            indices = np.argsort(self.sensitivity['s_tot_i'])[::-1]
            for i in indices:
                file.write(
                    '%s,%f,%f  \n' %
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

        # generate A matrix and evaluate the fitted surrogate
        a_matrix = self.calc_a(self.basis['multi-indices'])

        x = self.my_experiment.x_u
        data = np.dot(a_matrix, np.vstack((self.coefficients[:, 0])))

        # Concatenate along columns
        combined = np.hstack((x, data))

        # Prepare header
        header = self.my_experiment.my_data.stoch_data['names'] + [self.my_experiment.my_data.inputs['objective of interest']]

        # Open file and write everything
        output_file = os.path.join(
            self.my_experiment.my_data.path_res,
            "output_%s.csv" %
            self.my_experiment.my_data.inputs['objective of interest'])
        with open(output_file, 'w') as f:
            # Write header row
            f.write(','.join(header) + '\n')

            # Write data
            df = pd.DataFrame(combined)
            df.to_csv(f, header=False, index=False, lineterminator='\n')

        # generate the pdf
        density, bins = np.histogram(data, bins=100, density=1)
        centers = [(a + b) / 2 for a, b in zip(bins[::1], bins[1::1])]
        
        total = np.array([centers,density]).transpose()
        
        df1 = pd.DataFrame(total, columns=[self.my_experiment.my_data.inputs['objective of interest'],'probability density'])
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               ("data_pdf_%s.csv"
                                % (self.my_experiment.my_data.inputs[
                                    'objective of interest']))), "w") as f:
            df1.to_csv(f, index=False, lineterminator='\n')
                
        # generate the cdf
        cdf = np.cumsum(density * np.diff(bins))

        total = np.array([centers,cdf]).transpose()
        
        df1 = pd.DataFrame(total, columns=[self.my_experiment.my_data.inputs['objective of interest'],'cumulative probability'])
        with open(os.path.join(self.my_experiment.my_data.path_res,
                               ("data_cdf_%s.csv"
                                % (self.my_experiment.my_data.inputs[
                                    'objective of interest']))), "w") as f:
            df1.to_csv(f, index=False, lineterminator='\n')
