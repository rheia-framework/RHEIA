"""
The :py:mod:`determine_stoch_des_space` module contains a class to read the
design_space and stochastic_space files. In addition, functions are present
to load the case object and to evaluate the characterization of the input
dictionary.
"""

import os
import sys
import collections
import warnings

def custom_formatwarning(msg, *args, **kwargs):
    """
    Customize the format of the warning message:
    Show only the message.

    """
    # ignore everything except the message
    return str(msg) + '\n'



def check_dictionary(run_dict, uq_bool=False):
    """

    This function evaluates if the items in the input dictionary are
    properly characterized.

    Parameters
    ----------
    run_dict : dict
        The input dictionary.
    uq_bool : bool, optional
        Boolean that mentions if uncertainty quantification is considered.
        The default is False.

    """

    warnings.formatwarning = custom_formatwarning
    
    rob = False

    if not isinstance(run_dict, collections.Mapping):
        raise TypeError('The input dictionary should be a dictionary.')

    requirements = ['case',
                    'n jobs',
                    'results dir',
                    ]

    if 'n jobs' not in run_dict:
        run_dict['n jobs'] = 1

    for key in requirements:
        try:
            run_dict[key]
        except BaseException:
            raise KeyError(
                '"%s" is missing in the input dictionary.' %
                key)

    if not isinstance(run_dict['n jobs'], int):
        raise TypeError(
            'The value of the key "n jobs" should be a positive integer.')

    if not isinstance(run_dict['results dir'], str):
        raise TypeError(
            'The value of the key "results dir" should be a string.')

    if not uq_bool:
        requirements += ['objectives',
                         'stop',
                         'population size',
                         'x0',
                         'cx prob',
                         'mut prob',
                         'eta',
                         ]

        if 'x0' not in run_dict:
            run_dict['x0'] = ('AUTO', 'LHS')
            warnings.warn("If no previous results existed, the starting "
                          "population is generated using Latin Hypercube sampling.")
        if 'cx prob' not in run_dict:
            run_dict['cx prob'] = 0.9
            warnings.warn("The crossover probability is defined automatically "
                          "at 0.9")
        if 'mut prob' not in run_dict:
            run_dict['mut prob'] = 0.1
            warnings.warn("The mutation probability is defined automatically "
                          "at 0.1")
        if 'eta' not in run_dict:
            run_dict['eta'] = 0.2
            warnings.warn("The eta parameter is defined automatically at 0.2")

        for key in requirements[3:]:
            try:
                run_dict[key]
            except BaseException:
                raise KeyError(
                    '"%s" is missing in the input dictionary.' %
                    key)

        if not isinstance(run_dict['objectives'], collections.Mapping):
            raise TypeError(
                'The value of the key "objectives" should be a dictionary.')

        if not list(
                run_dict['objectives'].keys())[0] == 'DET' and not list(
                run_dict['objectives'].keys())[0] == 'ROB':
            raise ValueError(
                """Please select "DET" or "ROB" as a key
                for the "objectives" value.""")
        elif not isinstance(list(run_dict['objectives'].values())[0], tuple):
            raise TypeError(
                """The value of the key "DET" or "ROB" should be a tuple
                with the weights for the objectives
                (1 for maximization, -1 for minimization).""")
        elif (not all(isinstance(weight, int) for weight in
                      list(run_dict['objectives'].values())[0])):
            raise TypeError('The weights should be equal to 1 or -1.')
        elif (not all(abs(weight) == 1 for weight in
                      list(run_dict['objectives'].values())[0])):
            raise ValueError('The weights should be equal to 1 or -1.')

        if not isinstance(run_dict['stop'], tuple):
            raise TypeError(
                """The value of the key "stop" should be a tuple
                   with two elements.""")

        if not run_dict['stop'][0] == 'BUDGET':
            raise ValueError(
                """The first element in the tuple related to the key "stop"
                   should be equal to "BUDGET".""")

        if not isinstance(run_dict['stop'][1], int):
            raise TypeError(
                """The second element in the tuple related to the key "stop"
                   should be a positive integer.""")

        if not isinstance(run_dict['population size'], int):
            raise TypeError(
                """The value of the key "population number"
                   should be a positive integer.""")

        if not isinstance(run_dict['x0'], tuple):
            raise TypeError(
                """The value of the key "x0" should be a tuple
                   with two elements.""")

        if (not run_dict['x0'][0] == 'AUTO' and not
                run_dict['x0'][0] == 'CUSTOM'):
            raise ValueError(
                """The first element in the tuple related to the key "x0"
                   should be equal to "AUTO" or "CUSTOM".""")
        elif (run_dict['x0'][0] == 'AUTO' and not
              run_dict['x0'][1] == 'RANDOM' and not
              run_dict['x0'][1] == 'LHS'):
            raise ValueError(
                """When selecting "AUTO", the second element in the tuple
                   related to the key "x0" should be equal
                   to "RANDOM" or "LHS".""")

        if (not run_dict['x0'][0] == 'AUTO' and not
                run_dict['x0'][0] == 'CUSTOM'):
            raise ValueError(
                """The first element in the tuple related to the key "x0"
                   should be equal to "AUTO" or "CUSTOM".""")

        if not isinstance(run_dict['cx prob'], float):
            raise TypeError(
                """The value of the key "cx prob" should be a positive float,
                   lower or equal to 1.""")
        elif run_dict['cx prob'] > 1.:
            raise ValueError(
                """The value of the key "cx prob" should be a positive float,
                   lower or equal to 1.""")

        if not isinstance(run_dict['mut prob'], float):
            raise TypeError(
                """The value of the key "mut prob" should be a positive float,
                   lower or equal to 1.""")
        elif run_dict['mut prob'] > 1.:
            raise ValueError(
                """The value of the key "mut prob" should be a positive float,
                   lower or equal to 1.""")

        if not isinstance(run_dict['eta'], float):
            raise TypeError(
                """The value of the key "eta" should be a positive float,
                   lower or equal to 1.""")
        elif run_dict['eta'] > 1.:
            raise ValueError(
                """The value of the key "eta" should be a positive float,
                   lower or equal to 1.""")

        if list(run_dict['objectives'].keys())[0] == 'ROB':
            rob = True

    else:

        requirements += ['create only samples',
                         'objective of interest',
                         'draw pdf cdf',
                         ]

        if 'create only samples' not in run_dict:
            run_dict['create only samples'] = False

        if 'draw pdf cdf' not in run_dict:
            run_dict['draw pdf cdf'] = [False]

        for key in requirements[-3:]:
            try:
                run_dict[key]
            except BaseException:
                raise KeyError(
                    '"%s" is missing in the input dictionary.' %
                    key)

        if not isinstance(run_dict['create only samples'], bool):
            raise TypeError(
                """The value of the key "create only samples"
                should be a True or False.""")

        if not isinstance(run_dict['draw pdf cdf'], list):
            raise TypeError(
                'The value of the key "draw pdf cdf" should be a list.')

        if not isinstance(run_dict['draw pdf cdf'][0], bool):
            raise TypeError(
                """The value of the first element in the list "draw pdf cdf"
                should be a True or False.""")

        if run_dict['draw pdf cdf'][0] and not len(
                run_dict['draw pdf cdf']) == 2:
            raise ValueError(
                """When creating the distribution is desired,
                   provide the number of samples in the form of
                   a positive integer as the second list item""")

            if not isinstance(run_dict['draw pdf cdf'][1], int):
                raise TypeError(
                    """When creating the distribution is desired,
                    provide the number of samples in the form of
                    a positive integer as the second list item""")

        if not any(name == run_dict['objective of interest']
                   for name in run_dict['objective names']):
            raise ValueError(
                """The value of the key "objective of interest"
                   should be a string equal to any element
                   in the list with key "objective names".""")

    if rob or uq_bool:

        requirements += ['pol order',
                         'sampling method',
                         'objective names',
                         ]

        if 'sampling method' not in run_dict:
            run_dict['sampling method'] = 'SOBOL'
            warnings.warn("The Sobol' sequence is selected as "
                          "sampling method.")

        for key in requirements[-3:]:
            try:
                run_dict[key]
            except BaseException:
                raise KeyError(
                    '"%s" is missing in the input dictionary.' %
                    key)

        if not isinstance(run_dict['pol order'], int):
            raise TypeError(
                """The value of the key "pol order" should be a
                   positive integer.""")

        if (not run_dict['sampling method'] == 'RANDOM' and not
                run_dict['sampling method'] == 'SOBOL'):
            raise ValueError(
                """The value of the key "sampling method"
                   should be equal to "RANDOM" or "SOBOL".""")

        if not isinstance(run_dict['objective names'], list):
            raise TypeError(
                """The value of the key "objective names"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")

    if rob:
        requirements += ['objective of interest']

        try:
            run_dict[requirements[-1]]
        except BaseException:
            raise KeyError(
                '"%s" is missing in the input dictionary.' %
                requirements[-1])

        if not isinstance(run_dict['objective of interest'], list):
            raise TypeError(
                """The value of the key "objective of interest"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")

        if not isinstance(run_dict['objective names'], list):
            raise TypeError(
                """The value of the key "objective names"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")
        else:
            for name_qoi in run_dict['objective of interest']:
                if not any(name == name_qoi
                           for name in run_dict['objective names']):
                    raise ValueError(
                        """The value of the key "objective of interest"
                           should be a list with strings equal to any element
                           in the list with key "objective names".""")

        if not len(list(run_dict['objectives'].values())[
                   0]) == 2 * len(run_dict['objective of interest']):
            raise ValueError(
                """The number of objectives of interest should be equal to
                           two times the number of weigths. Note that the
                           weigths with uneven indices relate to the standard
                           deviation of the quantity of interest and should
                           therefore be equal to -1 (minimization).""")

    for elem in list(run_dict.keys()):
        if requirements.count(elem) == 0:
            raise ValueError("""" "%s" does not belong in the dictionary.
                                 Try renaming it, such that it corresponds to
                                 an expected dictionary item, or consider
                                 removing it. """ % elem)


def load_case(run_dict, design_space, uq_bool=False, create_only_samples=False):
    """
    For the selected case, the design variables and model parameters
    are loaded based on information from :file:`design_space`.
    In addition, the class object for the selected case is instantiated.

    Parameters
    ----------
    run_dict : dict
        The dictionary with information on the uncertainty quantification.
    design_space : string
        The design_space filename.
    uq_bool : bool, optional
        Indicates if uncertainty quantification is performed.
        The default is False.
    create_only_samples : bool, optional
        Indicates if only samples need to be created. The default is False.

    Returns
    -------
    space_obj : object
        The object with information on the design space and stochastic space
    eval_func : function
        The evaluate function that evaluates the system model.
    params : list, optional
        List with fixed data, used during model evaluation.
    """

    # define case directory
    case_path = os.path.join(os.path.split(os.path.dirname(
                             os.path.abspath(__file__)))[0],
                             'CASES')

    # enlist folders in cases directory
    dir_list = [item for item in os.listdir(case_path) if os.path.isdir(
                os.path.join(
                    case_path,
                    item))]

    # determine the evaluation type
    if uq_bool:
        eval_type = 'UQ'
    elif 'ROB' in run_dict['objectives'].keys():
        eval_type = 'ROB'
    elif 'DET' in run_dict['objectives'].keys():
        eval_type = 'DET'
    else:
        raise ValueError('Error in name of optimization type!')

    # get the StochasticDesignSpace object
    space_obj = StochasticDesignSpace(eval_type,
                                      run_dict['case'],
                                      design_space)

    # in case of optimization, attach objectives
    if not uq_bool:
        space_obj.attach_objectives(run_dict['objectives'])

    if not run_dict['case'] in dir_list:
        raise ValueError('Missing folder: %s folder is not found!' %
                         run_dict['case'])

    if not os.path.isfile(os.path.join(case_path, run_dict['case'],
                                       'case_description.py')):
        raise ValueError('Missing file: case_description.py not found!')

    if not create_only_samples:
        sys.path.insert(0, os.path.join(case_path, run_dict['case']))
        import case_description

        # determine the evaluate function from the considered case directory
        eval_func = case_description.evaluate

        # if fixed parameters for model evaluation are present,
        # get the params list
        if hasattr(case_description, 'set_params'):
            params = case_description.set_params()
        else:
            params = []

    else:
        eval_func = None
        params = []
        
    params.append({'results dir': run_dict['results dir']})

    return space_obj, eval_func, params


class StochasticDesignSpace(object):
    """
    Class which creates an object which characterizes the stochastic
    design space for the system model evaluation.

    Parameters
    ----------
    opt_type : string
       Defines the type run. Available types = ['DET', 'ROB', 'UQ']
    case : string
        The case name.
    design_space : string
        the design_space file name considered for evaluation

    """

    def __init__(self, opt_type, case, design_space):

        self.case = case
        self.design_space = design_space

        path = os.path.dirname(os.path.abspath(__file__))
        self.case_path = os.path.join(path, case)
        self.opt_type = opt_type
        self.l_b = []
        self.u_b = []
        self.par_dict = {}
        self.var_dict = {}
        self.upar_dict = {}
        self.obj = None

        self.read_design_space()
        self.read_stochastic_space()

    def read_design_space(self):
        """
        Reads the :file:`design_space.csv` file and extracts the deterministic
        values for the parameters and the bounds for the design variables.

        """

        path_to_read = os.path.join(self.case_path, self.design_space)

        # check if design_space file exists
        if not os.path.isfile(path_to_read):
            raise ValueError(
                """Missing file: "design_space.csv" or the name of the case does
                   not exist. Make sure that the name of the case is equal to
                   the name of the folder in CASES.""")

        with open(path_to_read, 'r') as file:
            for line in file:
                tmp = line.split(",")

                # store model parameter names and deterministic value
                if tmp[1] == 'par':
                    if len(tmp) > 3 and tmp[3] != '\n' and tmp[3] != '':
                        raise IndexError(
                            """ Wrong characterization of the parameter %s.
                                Is it supposed to be a variable?
                                Change "par" into "var".""" % tmp[0])

                    self.par_dict[tmp[0]] = float(tmp[2])

                # store design variable names and bounds
                elif tmp[1] == 'var':
                    if len(tmp) < 4:
                        raise IndexError(
                            """ Wrong characterization of the design variable %s.
                                Is it supposed to be a parameter?
                                Change "var" into "par".""" %
                            tmp[0])
                    self.var_dict[tmp[0]] = [float(tmp[2]), float(tmp[3])]
                    self.l_b.append(float(tmp[2]))
                    self.u_b.append(float(tmp[3]))
                    if float(tmp[2]) >= float(tmp[3]):
                        raise NameError(
                            """The lower bound is equal or greater than the
                            upper bound of %s""" % tmp[0])
                else:
                    raise ValueError(
                        """ The line %s does not mention if the
                            characterization corresponds to a parameter ("par")
                            or a design variable ("var").""" % tmp)
        self.n_dim = len(self.var_dict.keys())
        self.n_par = len(self.par_dict.keys())

        # check if the optimization case contains design variables
        if not self.var_dict and 'UQ' not in self.opt_type:
            raise NameError(
                'Variable list is empty! Case cannot be instatiated!')

    def read_stochastic_space(self):
        """
        Reads the :file:`stochastic_space` file and extracts the name,
        distribution and deviation from the mean of the stochastic parameters.

        """

        path_to_read = os.path.join(self.case_path, 'stochastic_space.csv')

        if any(
            x in self.opt_type for x in [
                'ROB',
                'UQ']) and os.path.isfile(path_to_read):

            # get characteristics of stochastic parameters
            with open(path_to_read, 'r') as file:
                for line in file:
                    tmp = line.split(",")
                    self.upar_dict[tmp[0]] = [tmp[1], tmp[2], float(tmp[3])]

        elif (any(x in self.opt_type for x in ['ROB', 'UQ'])
              and not os.path.isfile(path_to_read)):

            raise NameError(
                """Missing file: "stochastic_space.csv".
                   Uncertainty matrix cannot be built!""")

    def attach_objectives(self, obj):
        """
        Attaches the objectives of the configured optimization run.

        Parameters
        ----------
        obj : dict
            the objectives for the optimization run

        """

        self.obj = obj

    def convert_into_dictionary(self, x_i):
        """
        Convert the input sample for model evaluation into a dictionary.

        Parameters
        ----------
        x_i : ndarray
            Input sample for model evaluation.
            x_i = [xp1, xp2, ..., xpm, xd1, xd2, ..., xdn]
            Parameters are included for uncertainty quantification reasons.

        Returns
        -------
        {**parameters, **inputs} : dict
            The dictionary with the variable and parameter names as keys,
            and the input sample values as values.

        """

        # include parameters in first dictionary
        parameters = dict(zip(self.par_dict.keys(), x_i[:len(self.par_dict)]))
        if self.var_dict:
            # when design variables are present, include them in second dict
            inputs = dict(zip(self.var_dict.keys(), x_i[-len(self.var_dict):]))
        else:
            inputs = {}

        # concatenate two dict and return
        return {**parameters, **inputs}
