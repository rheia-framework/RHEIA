"""
The :py:mod:`genetic_algorithms` module contains a class that performs
the NSGA-II optimization algorithm.
"""

import multiprocessing as mp
import os
import random
from functools import partial
from copy import deepcopy
import numpy as np
from deap import creator, base, tools
import rheia.UQ.pce as uq


class NSGA2:
    """

    The NSGAII class includes methods to perform a deterministic and
    robust optimization using NSGA-II optimizer.

    Parameters
    ----------
    run_dict : dict
        The dictionary with user-defined values for
        the characterization of the optimization.
    space_obj : object
        The object with information on the design space
        and stochastic space
    start_from_last_gen : bool
      Boolean that indicates if results exist in the considered result
      directory.
    params : list
        Fixed parameters used during the model evaluation.

    """

    def __init__(self, run_dict, space_obj, start_from_last_gen, params):
        self.run_dict = run_dict
        self.space_obj = space_obj
        self.start_from_last_gen = start_from_last_gen
        self.params = params
        self.my_experiment = None
        self.objective_position = None
        self.opt_res_dir = os.path.join(os.path.abspath(
                                        os.path.join(
                                            os.path.dirname(__file__),
                                            '..')),
                                        'RESULTS',
                                        space_obj.case,
                                        list(run_dict['objectives'].keys())[0],
                                        run_dict['results dir'],
                                        )
    ########################
    # result files methods #
    ########################

    def init_opt(self):
        '''

        Initialization of the results directory and writing of the
        column names in the STATUS file.

        '''

        # check if results directory exists. If not, create one
        if not os.path.exists(self.opt_res_dir):
            os.makedirs(self.opt_res_dir)

        # initialize the STATUS file
        msg = '%8s%8s' % ('ITER', 'EVALS') + '\n' + '%8s%8s' % (0, 0)
        self.write_status(msg)

    def parse_status(self):
        """

        This method extracts the current generation number
        and the number of model evaluations performed.

        Returns
        -------
        ite : int
            The number of generations performed.
        evals : int
            The number of model evaluations performed.

        """

        # the directory of the STATUS file
        stat_dir = os.path.join(self.opt_res_dir, 'STATUS')

        with open(stat_dir, 'rb') as file:

            # read header
            file.readline()

            # read generation number and current number of evaluations
            line = file.readlines()[-1]

            ite, evals = int(line.split()[0]), int(line.split()[1])

        return ite, evals

    def write_status(self, msg):
        """

        A message is appended to the STATUS file.
        This message consists of the generation number
        and computational budget spent.

        Parameters
        ----------
        msg : str
            The message, to be appended to the STATUS file.

        """

        # the directory of the STATUS file
        stat_dir = os.path.join(self.opt_res_dir, 'STATUS')

        # add the message to the STATUS file
        file = open(stat_dir, 'a')
        file.write(msg + '\n')
        file.close()

    def append_points_to_file(self, individuals, filename):
        """

        This function is used to append the result (population or fitness)
        to the corresponding file.

        Parameters
        ----------
        individuals : list
            The values to be appended to the file.
        filename : str
            The filename where the points
            should be appended.

        """

        # the directory to the file in the results directory
        file_dir = os.path.join(self.opt_res_dir, filename)

        with open(file_dir, 'a') as file:

            for n_in in individuals:

                for item in n_in:

                    file.write('%.10f ' % item)

                file.write('\n')

            file.write('- \n')

    ##########################################
    # sample creation and evaluation methods #
    ##########################################

    def define_samples_to_eval(self, pop):
        """

        Defines the set of samples considered for evaluation.
        Adds the list of parameters to the initial samples
        drawn form the desing space.

        Parameters
        ----------
        pop : list
            The set of samples (population form, only variables)
            considered for evaluation

        Returns
        -------
        samples_to_eval : list
            The set of samples considered for evaluation
            samples_to_eval = [[p11, p12, ..., p1N, x11, x12, ..., x1M],
                              [p21, p22, ..., p2N, x21, x22, ..., x2M],
                              ...
                              [pk1, pk2, ..., pkN, xk1, xk2, ..., xkM],
                              ]

           where N is the number of parameters
                 M is the number of variables
                 k is the number of samples

        unc_samples_to_eval : list
            The set of samples considered for uncertainty quantification

        """

        # Attach parameters list to the given set of samples
        temp = np.tile(list(self.space_obj.par_dict.values()), (len(pop), 1))

        # check if samples in population have the appropriate length
        for pop_sample in pop:
            if len(pop_sample) != len(self.space_obj.var_dict):
                raise NameError(
                    """A sample in the list of starting design samples
                       does not match the number of design variables.""")

        # store the samples
        temp_samples = np.hstack((temp, np.array(pop))).tolist()
        if 'DET' in self.space_obj.opt_type:
            # return the sample, with deterministic values for the parameters
            # and values for the design variables chosen by the optimizer
            return temp_samples, []

        if 'ROB' in self.space_obj.opt_type:
            # in case of robust optimization, the values for the parameters
            # need to be defined by the PCE algorithm

            # list for the deterministic params and the stochastic params
            samples_to_eval, unc_samples_to_eval = [], []

            # for each design sample, generate the random values for the
            # stochastic params
            for sample in temp_samples:
                my_data = uq.Data(self.run_dict, self.space_obj)

                # acquire information on stochastic parameters
                my_data.read_stoch_parameters(
                    var_values=sample[-len(self.space_obj.var_dict):])

                # array for the number of quantities of interest considered
                self.objective_position = np.zeros(
                    len(self.run_dict['objective names']))
                for index, obj in enumerate(
                        self.run_dict['objective of interest']):

                    # for each quantity of interest, determine its position
                    # in the list
                    self.objective_position[index] = self.run_dict[
                        'objective names'].index(obj)

                # generate a random experiment for the quantity of interest
                self.my_experiment = uq.RandomExperiment(
                    my_data, self.objective_position)

                # create uniform/gaussian distributions and corresponding
                # orthogonal polynomials
                self.my_experiment.create_distributions()

                # determine the number of random samples needed for PCE
                self.my_experiment.n_terms()

                n_samples_unc = self.my_experiment.n_samples

                # no previous results available on the samples for this
                # specific PCE
                self.my_experiment.x_prev = np.array([])
                self.my_experiment.y_prev = np.array([])

                # create a design of experiment for the stochastic params
                self.my_experiment.create_samples(
                    size=self.my_experiment.n_samples)

                # the samples to construct the PCE
                unc_samples = self.my_experiment.x_u

                # produce the final set to evaluate
                temp_unc_samples = np.tile(sample, (n_samples_unc, 1))
                for j, elem in enumerate(self.space_obj.upar_dict.keys()):
                    # find the index of the parameter or variable in the list
                    par_var_list = list(self.space_obj.par_dict.keys())
                    par_var_list += list(self.space_obj.var_dict.keys())
                    unc_idx = par_var_list.index(elem)

                    temp_unc_samples[:, unc_idx] = self.my_experiment.x_u[:, j]

                # convert ndarray to list
                temp_unc_samples = list(temp_unc_samples)

                # add samples to the overall database
                samples_to_eval += temp_unc_samples
                unc_samples_to_eval += list(unc_samples)

            return samples_to_eval, unc_samples_to_eval

    def evaluate_samples(self, samples):
        """

        Evaluation of the set of samples. If the number of jobs
        is larger than 1, parallel processing of the samples
        is considered.

        Parameters
        ----------
        samples : list
            The set of samples considered for evaluation.

        Returns
        -------
        fitness : list
            The calculated fitness for the set of samples.

        """

        eval_dict = []

        # convert the list for each sample into a dictionary
        for sample in samples:
            eval_dict.append(self.space_obj.convert_into_dictionary(sample))

        # linear processing
        if self.run_dict['n jobs'] == 1:
            fitness = []
            for index, sample in enumerate(eval_dict):
                # evaluate the sample dictionary in the evaluate function
                # provide also the index of the sample in the list of samples
                fitness.append(self.run_dict['evaluate']((index, sample),
                                                         params=self.params))

        else:
            # multiprocessing of the sample evaluations
            # provide also the index of the sample in the list of samples
            pool = mp.Pool(processes=self.run_dict['n jobs'])
            fitness = pool.map(partial(self.run_dict['evaluate'],
                                       params=self.params),
                               enumerate(eval_dict))
            pool.close()

        return fitness

    def assign_fitness_to_population(self, pop, fitness, unc_samples):
        """

        Assigns the calulated fitness to the corresponding samples. In case
        of robust optimization, a PCE is constructed first, based on the
        random samples and the corresponding model output. Thereafter, the
        mean and standard deviation are extracted for the quantities of
        interest

        Parameters
        ----------
        pop : list
            Set of samples in population format, considered for evaluation.
        fitness : list
            The calculated fitness.
        unc_samples : list
            The set of samples for uncertainty quantification.

        Returns
        -------
        pop : list
            The population, appended with the fitness values.

        """

        if 'DET' in self.space_obj.opt_type:

            # in case of deterministic optimization, append the deterministic
            # result from the sample evaluations to the population
            for i, sample in enumerate(pop):
                pop[i].fitness.values = fitness[i]

        elif 'ROB' in self.space_obj.opt_type:
            # in case of robust optimization, the mean and standard deviation
            # need to be quantified based on the deterministic result from the
            # sample evaluations

            # Attach parameters list to the sample
            temp = np.tile(list(self.space_obj.par_dict.values()),
                           (len(pop), 1))
            for i in range(len(temp)):

                # store the fitness values temporarily
                temp_fitness = fitness[:self.my_experiment.n_samples]

                # remove these fitness values, as they correspond to a
                # deterministic evaluation of the samples
                del fitness[:self.my_experiment.n_samples]

                # store the random samples in the PCE instance variable
                # for previously generated samples
                self.my_experiment.x_prev = np.array(
                    unc_samples[i * self.my_experiment.n_samples:(i + 1) *
                                self.my_experiment.n_samples])

                # create a design of experiment based on these
                # previously generated samples
                self.my_experiment.create_samples()

                # for each quantity of interest, generate a PCE
                fitness_values = []
                for obj_position in self.objective_position:
                    if obj_position > len(temp_fitness[0]) - 1:
                        raise TypeError(""" The objective "%s" falls out
                                            of the range of predefined
                                            quantities of interest.
                                            Only %i outputs are returned
                                            from the model""" %
                                        (self.run_dict['objective names'][
                                            int(obj_position)],
                                         len(temp_fitness[0])))

                    # add the fitness values to the PCE instance variable
                    # that stores the deterministic results of the sample
                    # evaluations
                    self.my_experiment.y = np.array(temp_fitness)[
                        :, int(obj_position)].reshape(-1, 1)

                    # create a PCE object
                    my_pce = uq.PCE(self.my_experiment)

                    # evaluate the PCE
                    my_pce.run()

                    # Assign the mean and standard deviation as fitness values
                    fitness_values.append(my_pce.moments['mean'])
                    fitness_values.append(
                        np.sqrt(my_pce.moments['variance']))

                pop[i].fitness.values = fitness_values

        return pop

    def read_doe(self, doe_dir):
        """

        Read in the design of experiment.

        Parameters
        ----------
        doe_dir : str
            The directory of the file with the population.

        Returns
        -------
        doe : list
            The Design Of Experiments, i.e. the samples
            to be evaluated in the system model.

        """

        doe = []
        file = open(doe_dir, 'rb')

        # Read doe points
        for line in file:
            doe.append([float(i) for i in line.split()])

            # test if the doe points situate in between the variable bounds
            violate = False
            for index, elem in enumerate(doe[-1]):
                if (elem < self.space_obj.l_b[index] or
                        elem > self.space_obj.u_b[index]):
                    violate = True

            if violate:
                raise TypeError("""Design sample %s violates the
                                   design variable bounds. """ % str(doe[-1]))

        # test if the doe size equals the number of the population size
        # provided in the optimization dictionary
        if len(doe) != self.run_dict['population size']:
            raise TypeError(
                """The number of design samples in the starting population
                   file does not match with the population number provided
                   in the dictionary item "population number". """)

        return doe

    def eval_doe(self):
        """

        Evaluation of the Design Of Experiments (DoE). First, the DoE
        is read from the corresponding population file and stored in
        a list. Then, the design samples are evaluated in the model.
        When a population and fitness values are available from a previous
        optimization run, the method stores the final population and fitness
        values from this previous run. Finally, the population and
        corresponding fitness values are stored in the population and
        fitness file in the results directory, respectively.

        Returns
        -------
        current_pop : list
            The current population.

        """

        # DoE folder path
        path_doe = os.path.join(os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '..')),
            'OPT',
            'INPUTS',
            self.space_obj.case,
            '%iD' %
            self.space_obj.n_dim)

        # the DoE file
        file_doe = 'DOE_n%i' % self.run_dict['population size']

        # READ DoE file
        doe = self.read_doe(os.path.join(path_doe, file_doe))

        # evaluate DoE
        n_eval = 0
        current_pop = []
        for i, sol in enumerate(doe):

            # append the instance Individual
            current_pop.append(creator.Individual(sol))

        if not self.start_from_last_gen:
            # in case there is no information from previous results

            # create samples to evaluate
            individuals_to_eval, unc_samples = self.define_samples_to_eval(
                current_pop)

            # evaluate the samples
            fitnesses = self.evaluate_samples(individuals_to_eval)
            n_eval += len(individuals_to_eval)

            # assign fitness to the initial population
            current_pop = self.assign_fitness_to_population(current_pop,
                                                            fitnesses,
                                                            unc_samples)

            # update the population and fitness files
            self.append_points_to_file(current_pop, 'population')

            self.append_points_to_file([x.fitness.values for x in current_pop],
                                       'fitness')

            # update the STATUS file
            self.write_status('%8i%8i' % (1, n_eval))

        else:
            # get the population and fitness from the last available generation
            with open(os.path.join(self.opt_res_dir, 'fitness')) as file:

                # Read DOE points
                output = []
                for line in file.readlines()[-len(doe) - 1:-1]:
                    output.append([float(i) for i in line.split()])

            # add the fitness to the corresponding sample in the population
            for i, pop in enumerate(current_pop):
                pop.fitness.values = output[i]

        return current_pop

    ###################
    # NSGA-II methods #
    ###################

    def nsga2_1iter(self, current_pop):
        """

        Run one iteration of the NSGA-II optimizer.
        Based on the crossover and mutation probabilities, the
        offsprings are created and evaluated. Based on the fitness
        of these offsprings and of the current population,
        a new population is created.

        Parameters
        ----------
        current_pop : list
            The initial population.

        Returns
        -------
        new_pop : list
            The updated population.

        """

        # create offspring, clone of current population
        offspring = [deepcopy(ind) for ind in current_pop]

        # perform crossover
        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < self.run_dict['cx prob']:

                # apply crossover operator
                # in place editing of individuals
                tools.cxSimulatedBinaryBounded(ind1, ind2,
                                               self.run_dict['eta'],
                                               self.space_obj.l_b,
                                               self.space_obj.u_b)

                # set the fitness to an empty tuple of modified individuals
                del ind1.fitness.values
                del ind2.fitness.values

        # perform mutation
        for mutant in offspring:

            if random.random() < self.run_dict['mut prob']:

                # apply mutation operator
                # in place editing of individuals
                tools.mutPolynomialBounded(mutant, self.run_dict['eta'],
                                           self.space_obj.l_b,
                                           self.space_obj.u_b,
                                           self.run_dict['mut prob'])

                # set fitness to empty tuple of modified individuals
                del mutant.fitness.values

        # evaluate the modified individuals
        n_eval = 0
        invalid_indices = []
        invalid_fit_individuals = []
        for i, ind in enumerate(offspring):

            if not ind.fitness.valid:

                invalid_indices.append(i)
                invalid_fit_individuals.append(ind)

        if len(invalid_fit_individuals) > 0:

            # create samples to evaluate
            individuals_to_eval, unc_samples = self.define_samples_to_eval(
                invalid_fit_individuals)

            # evaluate samples
            fitnesses = self.evaluate_samples(individuals_to_eval)
            n_eval += len(individuals_to_eval)

            # assign fitness to the orginal samples list
            individuals_to_assign = self.assign_fitness_to_population(
                invalid_fit_individuals,
                fitnesses,
                unc_samples)

            # construct offspring list
            for i, ind in zip(invalid_indices, individuals_to_assign):

                offspring[i] = deepcopy(ind)

        # select next population using the NSGA-II operator
        new_pop = tools.selNSGA2(current_pop + offspring, len(current_pop))

        # update the population and fitness files
        self.append_points_to_file(new_pop, 'population')

        fitness_values = [x.fitness.values for x in new_pop]
        self.append_points_to_file(fitness_values, 'fitness')

        # update the STATUS file
        ite, evals = self.parse_status()
        self.write_status('%8i%8i' % (ite + 1, evals + n_eval))

        return new_pop

    def run_optimizer(self):
        """

        Run an optimization using the NSGA-II algorithm.

        """

        weigth = list(self.space_obj.obj.values())[0]
        creator.create('Fitness', base.Fitness, weights=weigth)
        creator.create('Individual', list, fitness=creator.Fitness)

        if not self.start_from_last_gen:
            self.init_opt()

        # set up the STATUS file
        ite, init_evals = self.parse_status()

        # set up the logbook
        logbook = tools.Logbook()
        logbook.header = 'gen', 'evals', 'min', 'max'
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register('min', np.min, axis=0)
        stats.register('max', np.max, axis=0)

        # evaluate the DoE
        temp_output_pop = self.eval_doe()

        # run the generations
        evals = 0
        while evals < self.run_dict['stop'][1]:

            # perform one NSGA-II iteration
            temp_output_pop = self.nsga2_1iter(temp_output_pop)

            # update the evaluations counter
            ite, current_evals = self.parse_status()
            evals = current_evals - init_evals

            # update the logbook record
            record = stats.compile(temp_output_pop)
            logbook.record(gen=ite, evals=evals, **record)
            print(logbook.stream)


def return_opt_methods():
    """

    Returns the names of the available genetic optimizers.

    Returns
    -------
    name_list : list
        A list with the names of the available genetic optimizers.

    """

    name_list = ['NSGA2']

    return name_list


def return_opt_obj(name):
    """

    Returns the optimizer fuction object

    Parameters
    ----------
    name : str
        Name of the optimizer.

    Returns
    -------
    object
        The optimizer function object.

    """

    switcher = {'NSGA2': NSGA2,
                }

    return switcher.get(name)
