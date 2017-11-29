acf
--------

Installation:

Clone this repository and cd into it.
    >>> python setup.py install

To measure a rotation period from a light curve:

    >>> from simple_acf import simple_acf
    >>> period, acf, lags, rvar, peaks = simple_acf(time, flux)

or::

    >>> from Kepler_ACF import corr_run
    >>> period, period_error = corr_run(time, flux, flux_err, id, savedir,
        saveplot=True)
