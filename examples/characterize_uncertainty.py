"""
Characterize uncertain parameters from two percentiles (Cook's method).

Uncertainty is often available as two percentiles rather than as a mean and
a deviation directly, e.g. a declared value at the 50th percentile and a
design value at the 99th percentile from a standard or specification. Two
percentiles pin down a two-parameter distribution exactly, so they can be
converted straight into the mean and deviation RHEIA's design_space.csv and
stochastic_space.csv expect.

To run this example:

1. Create a case folder under CASES (e.g. CASES/MY_CASE).
2. Copy CASES/NO_MODEL/uncertainty_characterization.csv into it and fill in
   one row per uncertain parameter:

       name,distribution,x_low,p_low,x_high,p_high
       th_cond,Lognormal,500,0.5,3000,0.99

   where x_low is the value at percentile p_low, x_high the value at
   percentile p_high (both percentiles as fractions between 0 and 1), and
   distribution is Uniform, Gaussian or Lognormal.
3. Set CASE below to the folder name and run this script.

Each row is fitted and written into that case's design_space.csv (as a
model parameter, holding the fitted mean) and stochastic_space.csv (holding
the fitted deviation, with an absolute relation). Both files are created
when they do not exist yet; parameters already present are replaced, and
rows the characterization file does not mention are left untouched.

Picking good percentiles is usually iterative: start with WRITE = False and
SHOW = True, re-run as you tweak the values in the characterization file and
look at the plotted fit, then set WRITE = True once the fit looks right.
Nothing is written to disk at all while WRITE is False.

The example below reproduces the motivating use case: a material's thermal
conductivity, declared as 500 W/mK at the 50th percentile and 3000 W/mK at
the 99th percentile, fitted to a Lognormal so the value stays positive.
"""

import rheia.UQ.characterize_uncertainty as rheia_char

# the folder name under CASES, holding the
# uncertainty_characterization.csv file
CASE = 'MY_CASE'

# False: just fit/plot, nothing touches disk -- iterate on the characterization
# file until the fit looks right.
# True: actually write design_space.csv and stochastic_space.csv.
WRITE = False

# replace a parameter that is already characterized in the case;
# False raises on it instead
OVERWRITE = True

# render the fitted PDF with the two percentiles marked
PLOT = True

# display it interactively (blocks until the plot window is closed); also saved
# as '<name>_fit.png' in the case folder once WRITE is True
SHOW = True


if __name__ == '__main__':
    fitted = rheia_char.characterize_from_csv(
        CASE, overwrite=OVERWRITE, write=WRITE, plot=PLOT, show=SHOW)

    for name, (value, deviation) in fitted.items():
        print('%s: value=%.6g, deviation=%.6g' % (name, value, deviation))
