"""
The :py:mod:`determine_stoch_des_space` module contains a class to read the
design_space and stochastic_space files. In addition, functions are present
to load the case object and to evaluate the characterization of the input
dictionary.
"""

import os
import importlib.util
import sys
import collections.abc
import warnings
import pandas as pd

DESIGN_SPACE_COLUMNS = ['name', 'type', 'value', 'upper_bound']
STOCHASTIC_SPACE_COLUMNS = ['name', 'relation', 'distribution', 'deviation']


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

    if not isinstance(run_dict, collections.abc.Mapping):
        raise TypeError('The input dictionary should be a dictionary.')

    required = {'case', 'n jobs', 'results dir'}
    allowed = set(required)
    defaults = {'n jobs': 1}
    rob = False
    weights = None

    if uq_bool:
        required.update({'objective of interest', 'draw pdf cdf'})
        allowed.update({'create only samples', 'objective of interest',
                        'draw pdf cdf'})
        defaults.update({'create only samples': False,
                         'draw pdf cdf': [False]})
    else:
        required.update({'objectives', 'stop', 'population size'})
        allowed.update({'objectives', 'stop', 'population size', 'x0',
                        'cx prob', 'mut prob', 'eta'})
        defaults.update({'x0': ('AUTO', 'LHS'),
                         'cx prob': 0.9,
                         'mut prob': 0.1,
                         'eta': 0.2})

    for key, value in defaults.items():
        if key in run_dict:
            continue

        run_dict[key] = value
        if key == 'x0':
            warnings.warn("If no previous results existed, the starting "
                          "population is generated using Latin Hypercube sampling.")
        elif key == 'cx prob':
            warnings.warn("The crossover probability is defined automatically "
                          "at 0.9")
        elif key == 'mut prob':
            warnings.warn("The mutation probability is defined automatically "
                          "at 0.1")
        elif key == 'eta':
            warnings.warn("The eta parameter is defined automatically at 0.2")

    for key in required:
        if key not in run_dict:
            raise KeyError('"%s" is missing in the input dictionary.' % key)

    if (not isinstance(run_dict['n jobs'], int) or
            isinstance(run_dict['n jobs'], bool)):
        raise TypeError(
            'The value of the key "n jobs" should be a positive integer.')
    if run_dict['n jobs'] <= 0:
        raise ValueError(
            'The value of the key "n jobs" should be a positive integer.')

    if not isinstance(run_dict['results dir'], str):
        raise TypeError(
            'The value of the key "results dir" should be a string.')
    if run_dict['results dir'] == '':
        raise ValueError(
            'The value of the key "results dir" should not be empty.')

    if not uq_bool:
        if not isinstance(run_dict['objectives'], collections.abc.Mapping):
            raise TypeError(
                'The value of the key "objectives" should be a dictionary.')
        if len(run_dict['objectives']) != 1:
            raise ValueError(
                'The value of the key "objectives" should contain one run type.')

        opt_type = next(iter(run_dict['objectives'].keys()))
        weights = run_dict['objectives'][opt_type]
        rob = opt_type == 'ROB'

        if opt_type not in ['DET', 'ROB']:
            raise ValueError(
                """Please select "DET" or "ROB" as a key
                for the "objectives" value.""")
        if not isinstance(weights, tuple):
            raise TypeError(
                """The value of the key "DET" or "ROB" should be a tuple
                with the weights for the objectives
                (1 for maximization, -1 for minimization).""")
        if not weights:
            raise ValueError('The objective weights tuple should not be empty.')
        if not all(isinstance(weight, int) for weight in weights):
            raise TypeError('The weights should be equal to 1 or -1.')
        if not all(abs(weight) == 1 for weight in weights):
            raise ValueError('The weights should be equal to 1 or -1.')

        if not isinstance(run_dict['stop'], tuple):
            raise TypeError(
                """The value of the key "stop" should be a tuple
                   with two elements.""")
        if len(run_dict['stop']) != 2:
            raise ValueError(
                """The value of the key "stop" should be a tuple
                   with two elements.""")
        if run_dict['stop'][0] != 'BUDGET':
            raise ValueError(
                """The first element in the tuple related to the key "stop"
                   should be equal to "BUDGET".""")
        if (not isinstance(run_dict['stop'][1], int) or
                isinstance(run_dict['stop'][1], bool)):
            raise TypeError(
                """The second element in the tuple related to the key "stop"
                   should be a positive integer.""")
        if run_dict['stop'][1] <= 0:
            raise ValueError(
                """The second element in the tuple related to the key "stop"
                   should be a positive integer.""")

        if (not isinstance(run_dict['population size'], int) or
                isinstance(run_dict['population size'], bool)):
            raise TypeError(
                """The value of the key "population number"
                   should be a positive integer.""")
        if run_dict['population size'] <= 0:
            raise ValueError(
                """The value of the key "population number"
                   should be a positive integer.""")

        if not isinstance(run_dict['x0'], tuple):
            raise TypeError(
                """The value of the key "x0" should be a tuple
                   with two elements.""")
        if len(run_dict['x0']) != 2:
            raise ValueError(
                """The value of the key "x0" should be a tuple
                   with two elements.""")
        if run_dict['x0'][0] not in ['AUTO', 'CUSTOM']:
            raise ValueError(
                """The first element in the tuple related to the key "x0"
                   should be equal to "AUTO" or "CUSTOM".""")
        if (run_dict['x0'][0] == 'AUTO' and
                run_dict['x0'][1] not in ['RANDOM', 'LHS']):
            raise ValueError(
                """When selecting "AUTO", the second element in the tuple
                   related to the key "x0" should be equal
                   to "RANDOM" or "LHS".""")

        for key in ['cx prob', 'mut prob', 'eta']:
            if not isinstance(run_dict[key], float):
                raise TypeError(
                    """The value of the key "%s" should be a positive float,
                       lower or equal to 1.""" % key)
            if run_dict[key] <= 0. or run_dict[key] > 1.:
                raise ValueError(
                    """The value of the key "%s" should be a positive float,
                       lower or equal to 1.""" % key)

    else:
        if not isinstance(run_dict['create only samples'], bool):
            raise TypeError(
                """The value of the key "create only samples"
                should be a True or False.""")

        if not isinstance(run_dict['draw pdf cdf'], list):
            raise TypeError(
                'The value of the key "draw pdf cdf" should be a list.')
        if not run_dict['draw pdf cdf']:
            raise ValueError(
                'The value of the key "draw pdf cdf" should not be empty.')
        if not isinstance(run_dict['draw pdf cdf'][0], bool):
            raise TypeError(
                """The value of the first element in the list "draw pdf cdf"
                should be a True or False.""")
        if run_dict['draw pdf cdf'][0] and len(run_dict['draw pdf cdf']) != 2:
            raise ValueError(
                """When creating the distribution is desired,
                   provide the number of samples in the form of
                   a positive integer as the second list item""")
        if run_dict['draw pdf cdf'][0]:
            n_draw_samples = run_dict['draw pdf cdf'][1]
            if isinstance(n_draw_samples, bool):
                raise TypeError(
                    """When creating the distribution is desired,
                    provide the number of samples in the form of
                    a positive integer as the second list item""")
            if isinstance(n_draw_samples, float) and n_draw_samples.is_integer():
                run_dict['draw pdf cdf'][1] = int(n_draw_samples)
                n_draw_samples = run_dict['draw pdf cdf'][1]
            if not isinstance(n_draw_samples, int):
                raise TypeError(
                    """When creating the distribution is desired,
                    provide the number of samples in the form of
                    a positive integer as the second list item""")
            if n_draw_samples <= 0:
                raise ValueError(
                    """When creating the distribution is desired,
                    provide the number of samples in the form of
                    a positive integer as the second list item""")

    if rob or uq_bool:
        required.update({'pol order', 'sampling method', 'objective names',
                         'uq method'})
        allowed.update({'pol order', 'sampling method', 'objective names',
                        'uq method', 'n samples'})

        if 'sampling method' not in run_dict:
            run_dict['sampling method'] = 'SOBOL'
            warnings.warn("The Sobol' sequence is selected as "
                          "sampling method.")
        if 'uq method' not in run_dict:
            run_dict['uq method'] = 'full'
            warnings.warn("A full conventional PCE is selected.")
        if run_dict['uq method'] == 'sparse' and 'n samples' not in run_dict:
            run_dict['n samples'] = 10
            warnings.warn("The number of sparse PCE samples is defined "
                          "automatically at 10.")

        for key in required:
            if key not in run_dict:
                raise KeyError(
                    '"%s" is missing in the input dictionary.' % key)

        if (not isinstance(run_dict['pol order'], int) or
                isinstance(run_dict['pol order'], bool)):
            raise TypeError(
                """The value of the key "pol order" should be a
                   positive integer.""")
        if run_dict['pol order'] <= 0:
            raise ValueError(
                """The value of the key "pol order" should be a
                   positive integer.""")

        if run_dict['uq method'] not in ['full', 'sparse']:
            raise ValueError(
                """The value of the key "uq method" should be equal to
                   "full" or "sparse".""")
        if run_dict['uq method'] == 'sparse':
            if (not isinstance(run_dict['n samples'], int) or
                    isinstance(run_dict['n samples'], bool)):
                raise TypeError(
                    """The value of the key "n samples" should be a
                       positive integer for sparse PCE.""")
            if run_dict['n samples'] <= 0:
                raise ValueError(
                    """The value of the key "n samples" should be a
                       positive integer for sparse PCE.""")

        if run_dict['sampling method'] not in ['RANDOM', 'SOBOL']:
            raise ValueError(
                """The value of the key "sampling method"
                   should be equal to "RANDOM" or "SOBOL".""")

        if not isinstance(run_dict['objective names'], list):
            raise TypeError(
                """The value of the key "objective names"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")
        if not run_dict['objective names']:
            raise ValueError(
                """The value of the key "objective names"
                   should not be empty.""")
        if not all(isinstance(name, str) for name in run_dict['objective names']):
            raise TypeError(
                """The value of the key "objective names"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")

    if uq_bool:
        if not isinstance(run_dict['objective of interest'], str):
            raise TypeError(
                """The value of the key "objective of interest"
                   should be a string.""")
        if run_dict['objective of interest'] not in run_dict['objective names']:
            raise ValueError(
                """The value of the key "objective of interest"
                   should be a string equal to any element
                   in the list with key "objective names".""")

    if rob:
        required.add('objective of interest')
        allowed.add('objective of interest')

        if 'objective of interest' not in run_dict:
            raise KeyError(
                '"objective of interest" is missing in the input dictionary.')
        if not isinstance(run_dict['objective of interest'], list):
            raise TypeError(
                """The value of the key "objective of interest"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")
        if not all(isinstance(name, str)
                   for name in run_dict['objective of interest']):
            raise TypeError(
                """The value of the key "objective of interest"
                   should be a list with the names of the
                   quantities of interest as elements in string format.""")
        for name_qoi in run_dict['objective of interest']:
            if name_qoi not in run_dict['objective names']:
                raise ValueError(
                    """The value of the key "objective of interest"
                       should be a list with strings equal to any element
                       in the list with key "objective names".""")

        if len(weights) != 2 * len(run_dict['objective of interest']):
            raise ValueError(
                """The number of objectives of interest should be equal to
                           two times the number of weigths. Note that the
                           weigths with uneven indices relate to the standard
                           deviation of the quantity of interest and should
                           therefore be equal to -1 (minimization).""")

    for elem in list(run_dict.keys()):
        if elem not in allowed:
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
        case_description_path = os.path.join(
            case_path, run_dict['case'], 'case_description.py')
        module_name = 'rheia.CASES.%s.case_description' % run_dict['case']
        spec = importlib.util.spec_from_file_location(
            module_name, case_description_path)
        case_description = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = case_description
        spec.loader.exec_module(case_description)

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
        params = None

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

        data = pd.read_csv(path_to_read, dtype=str, keep_default_na=False,
                           na_filter=False)
        missing_columns = [column for column in DESIGN_SPACE_COLUMNS
                           if column not in data.columns]
        if missing_columns:
            raise ValueError(
                """The design_space file should contain the columns: %s."""
                % ', '.join(DESIGN_SPACE_COLUMNS))

        data = data[DESIGN_SPACE_COLUMNS]
        data = data.apply(lambda col: col.str.strip())

        for row in data.itertuples(index=False):
            if row.name == '' or row.type == '' or row.value == '':
                raise ValueError(
                    """The design_space row %s is incomplete.""" %
                    (row.Index if hasattr(row, 'Index') else row,))

            if row.type == 'par':
                if row.upper_bound != '':
                    raise IndexError(
                        """ Wrong characterization of the parameter %s.
                            Is it supposed to be a variable?
                            Change "par" into "var".""" % row.name)

                self.par_dict[row.name] = float(row.value)

            elif row.type == 'var':
                if row.upper_bound == '':
                    raise IndexError(
                        """ Wrong characterization of the design variable %s.
                            Is it supposed to be a parameter?
                            Change "var" into "par".""" % row.name)

                lower = float(row.value)
                upper = float(row.upper_bound)
                self.var_dict[row.name] = [lower, upper]
                self.l_b.append(lower)
                self.u_b.append(upper)
                if lower >= upper:
                    raise NameError(
                        """The lower bound is equal or greater than the
                        upper bound of %s""" % row.name)
            else:
                raise ValueError(
                    """ The line %s does not mention if the
                        characterization corresponds to a parameter ("par")
                        or a design variable ("var").""" % (row,))
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

            data = pd.read_csv(path_to_read, dtype=str,
                               keep_default_na=False, na_filter=False)
            missing_columns = [column for column in STOCHASTIC_SPACE_COLUMNS
                               if column not in data.columns]
            if missing_columns:
                raise ValueError(
                    """The stochastic_space file should contain the columns: %s."""
                    % ', '.join(STOCHASTIC_SPACE_COLUMNS))

            data = data[STOCHASTIC_SPACE_COLUMNS]
            data = data.apply(lambda col: col.str.strip())

            for row in data.itertuples(index=False):
                if any(value == '' for value in row):
                    raise ValueError(
                        """The stochastic_space row %s is incomplete."""
                        % (row,))

                relation = row.relation
                distribution = row.distribution.capitalize()
                if relation not in ['absolute', 'relative']:
                    raise ValueError(
                        """The uncertainty relation for %s should be
                           "absolute" or "relative".""" % row.name)
                if distribution not in ['Uniform', 'Gaussian']:
                    raise ValueError(
                        """The uncertainty distribution for %s should be
                           "Uniform" or "Gaussian".""" % row.name)

                self.upar_dict[row.name] = [
                    relation, distribution, float(row.deviation)]

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
