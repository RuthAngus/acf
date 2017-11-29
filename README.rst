acf
--------

Installation:

Clone this repository and cd into it.
    >>> python setup.py install

To measure an ACF rotation period from a light curve:

    >>> from simple_acf import simple_acf
    >>> period, acf, lags, rvar, peaks = simple_acf(time, flux)

or::

    >>> from Kepler_ACF import corr_run
    >>> period, period_error = corr_run(time, flux, flux_err, id, savedir,
        saveplot=True)

To measure a Lomb-Scargle periodogram period::

    >>> from rotation import prot as prt
    >>> prot = prt(kepid, x, y, yerr, LC_DIR)
    >>> pgram_period, pgram_period_err = prot.pgram_ps(plot=True)

Where x, y, yerr are the time, flux and flux_err arrays.
LC_DIR is the path to the kplr downloaded light curves.
x, y, yerr and LC_DIR are optional, if you do not provide them, kplr will
download the light curve.

If plot=True in pgram_ps, dianostic plots will be made.
