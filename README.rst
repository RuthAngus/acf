acf
--------

To infer an age from a rotation period and (B-V) colour::

    >>> from simple_acf import simple_acf
    >>> period, acf, lags, rvar, peaks = simple_acf(time, flux)

or::

    >>> from Kepler_ACF import corr_run
    >>> period, period_error = corr_run(time, flux, flux_err, id, savedir,
        saveplot=True)
