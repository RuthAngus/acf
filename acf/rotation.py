"""
A system for measuring rotation periods.
This contains functions for measuring rotation.
"""

import os
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import filtering as flt

import kplr
client = kplr.API()

import kepler_data as kd
from astropy.stats import LombScargle


class prot(object):
    """
    Given a star object with a kepid or x, y and yerr values, measure the
    rotation period.
    __init__ downloads the light curve if it doesn't already exist.
    pgram_ps measures a periodogram rotation period.
    """

    def __init__(self, kepid=None, x=None, y=None, yerr=None,
                 LC_DIR="/Users/ruthangus/.kplr/data/lightcurves"):
        """
        params:
        ------
        kepid: (int)
            The KIC id.
        x: (array)
            The time array.
        y: (array)
            The flux array.
        yerr: (array)
            The flux uncertainty array.
        """
        self.kepid = kepid

        # If x, y and yerr are not provided, load them.
        if not np.array([x, y, yerr]).any():
            lc_path = os.path.join(LC_DIR, str(kepid).zfill(9))

            # If you don't have the light curves, download them.
            if not os.path.exists(lc_path):
                print("Downloading light curve...")
                star = client.star(kepid)
                star.get_light_curves(fetch=True, short_cadence=False)

            print("Loading light curve...")
            self.x, self.y, self.yerr = kd.load_kepler_data(lc_path)
        else:
            self.x, self.y, self.yerr = x, y, yerr

    def pgram_ps(self, filter_period=35., plot=False, clobber=False):
        """
        Measure a periodogram rotation period
        parameters:
        ----------
        filter_period: (float or None)
            if None, the lc is not filtered.
            if float, the lc is high-pass filtered with a cutoff of
            filter_period.
        plot: (bool)
            If true the periodogram is plotted and saved.
        clobber: (bool)
            If true any existing periodogram file is overwritten.

        returns:
        -------
        pgram_period: (float)
            The rotation period
        pgram_period_err: (float)
            The formal uncertainty on the period.

        Adds self.pgram_period, self.pgram_period.err, self.pgram and self.ps.
        """

        pgram_fname = "pgrams/{}_pgram".format(self.kepid)
        if not os.path.exists("pgrams"):
            os.mkdir("pgrams")
        if clobber:
            freq = np.linspace(1./100, 1./.1, 100000)
            ps = 1./freq

            if filter_period:
                filter_period = filter_period  # days
                fs = 1./(self.x[1] - self.x[0])
                lowcut = 1./filter_period
                yfilt = flt.butter_bandpass_filter(self.y, lowcut, fs, order=3,
                                                plot=False)

            else:
                yfilt = self.y*1

            print("Calculating periodogram")
            pgram = LombScargle(self.x, yfilt, self.yerr).power(freq)

            peaks = np.array([i for i in range(1, len(ps)-1) if pgram[i-1] <
                                pgram[i] and pgram[i+1] < pgram[i]])

            presults = pd.DataFrame({"periods": ps, "power": pgram})
            presults.to_csv("{}.csv".format(pgram_fname))
            print("saving pgram to {}.csv".format(pgram_fname))

        else:  # If not clobber, look for old result.
            print("looking for {}.csv".format(pgram_fname))
            if os.path.exists("{}.csv".format(pgram_fname)):  # If pgram already exists
                print("{}.csv found, loading pre-exisiting periodogram"
                      .format(pgram_fname))
                pr = pd.read_csv("{}.csv".format(pgram_fname))
                ps, pgram = pr.periods.values, pr.power.values

            else:  # If pgram does not exist.
                freq = np.linspace(1./100, 1./.1, 100000)
                ps = 1./freq

                filter_period = 35.  # days
                fs = 1./(self.x[1] - self.x[0])
                lowcut = 1./filter_period
                yfilt = flt.butter_bandpass_filter(self.y, lowcut, fs, order=3,
                                                plot=False)

                print("Calculating periodogram")
                pgram = LombScargle(self.x, yfilt, self.yerr).power(freq)

                presults = pd.DataFrame({"periods": ps, "power": pgram})
                presults.to_csv("{}.csv".format(pgram_fname))
                print("saving pgram to {}.csv".format(pgram_fname))

        peaks = np.array([i for i in range(1, len(ps)-1) if pgram[i-1] <
                            pgram[i] and pgram[i+1] < pgram[i]])

        pgram_period = ps[pgram == max(pgram[peaks])][0]
        if plot:
            pl.clf()
            pl.subplot(2, 1, 1)
            pl.plot(self.x-self.x[0], self.y, "k.", ms=3)
            pl.xlim(0, 50)
            pl.title("Period = {0:.2f} days".format(pgram_period))
            pl.subplot(2, 1, 2)
            pl.plot(ps, pgram)
            pl.axvline(pgram_period, color="orange", ls="--")
            pl.xlabel("Period (days)")
            pl.savefig(pgram_fname)
            print("saving plot as {}.png".format(pgram_fname))

        # Calculate the uncertainty.
        _freq = 1./pgram_period
        pgram_freq_err = self.calc_pgram_uncertainty(_freq)
        frac_err = pgram_freq_err/_freq
        pgram_period_err = pgram_period * frac_err

        self.pgram_period = pgram_period
        self.pgram_period_err = pgram_period_err
        self.pgram = pgram
        self.ps = ps
        return pgram_period, pgram_period_err

    def calc_pgram_uncertainty(self, freq):
        """
        Calculate the formal uncertainty on the periodogram period (1/freq).
        """
        phase, A = self.calc_phase_and_amp(freq)
        y_noise = self.y - A**2*np.sin(2*np.pi*freq*self.x + phase)
        sigma_n = np.var(y_noise)
        N, T = len(self.x), self.x[-1] - self.x[0]
        return 3 * np.pi * sigma_n / (2 * N**.5 * T * A)

    def calc_phase_and_amp(self, f):
        """
        Phase and amplitude calculation for the calc_pgram_uncertainty
        function.
        """
        AT = np.vstack((self.x, np.ones((3, len(self.y)))))
        ATA = np.dot(AT, AT.T)
        arg = 2*np.pi*f*self.x
        AT[-2, :] = np.sin(arg)
        AT[-1, :] = np.cos(arg)
        v = np.dot(AT[:-2, :], AT[-2:, :].T)
        ATA[:-2, -2:] = v
        ATA[-2:, :-2] = v.T
        ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)
        w = np.linalg.solve(ATA, np.dot(AT, self.y))
        A, B = w[-1], w[-2]
        phase = np.arctan(A/B)
        Amp = (np.abs(A) + np.abs(B))**.5
        return phase, Amp
