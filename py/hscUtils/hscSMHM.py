"""Parameterized models of the stellar mass - halo mass relation (SMHM)."""

from __future__ import division, print_function
from __future__ import absolute_import, unicode_literals

import numpy as np

__all__ = ['Leauthaud12']


class Leauthaud12():
    """
    Class for stellar mass halo mass relation from Leauthaud 2012.
    """

    def __init__(self, redshift=0.3, sigmod=1):
        """Initiallize the model."""
        if (redshift > 0.22) and (redshift <= 0.48):

            if sigmod == 1:

                self.param = {'logM1': 12.520, 'logM1Err': 0.037,
                              'logMs0': 10.916, 'logMs0Err': 0.020,
                              'beta': 0.457, 'betaErr': 0.009,
                              'delta': 0.566, 'deltaErr': 0.086,
                              'gamma': 1.530, 'gammaErr': 0.180,
                              'sigMs': 0.206, 'sigMsErr': 0.031,
                              'bCut': 1.47, 'bCutErr': 0.73,
                              'bSat': 10.62, 'bSatErr': 0.87,
                              'betaCut': -0.13, 'betaCutErr': 0.28,
                              'betaSat': 0.859, 'betaSatErr': 0.038}

            elif sigmod == 2:

                self.param = {'logM1': 12.518, 'logM1Err': 0.038,
                              'logMs0': 10.917, 'logMs0Err': 0.020,
                              'beta': 0.456, 'betaErr': 0.009,
                              'delta': 0.582, 'deltaErr': 0.083,
                              'gamma': 1.480, 'gammaErr': 0.170,
                              'sigMs': 0.192, 'sigMsErr': 0.031,
                              'bCut': 1.52, 'bCutErr': 0.79,
                              'bSat': 10.69, 'bSatErr': 0.89,
                              'betaCut': -0.11, 'betaCutErr': 0.29,
                              'betaSat': 0.860, 'betaSatErr': 0.039}

            else:
                raise KeyError("Wrong SIG_MOD choice!!")

        elif (redshift > 0.48) and (redshift <= 0.74):

            if sigmod == 1:

                self.param = {'logM1': 12.725, 'logM1Err': 0.032,
                              'logMs0': 11.038, 'logMs0Err': 0.019,
                              'beta': 0.466, 'betaErr': 0.009,
                              'delta': 0.610, 'deltaErr': 0.130,
                              'gamma': 1.950, 'gammaErr': 0.250,
                              'sigMs': 0.249, 'sigMsErr': 0.019,
                              'bCut': 1.65, 'bCutErr': 0.65,
                              'bSat': 9.04, 'bSatErr': 0.81,
                              'betaCut': 0.590, 'betaCutErr': 0.280,
                              'betaSat': 0.740, 'betaSatErr': 0.059}

            else:
                raise KeyError("Wrong SIG_MOD choice!!")

        elif (redshift > 0.74) and (redshift <= 1.00):

            if sigmod == 1:

                self.param = {'logM1': 12.722, 'logM1Err': 0.027,
                              'logMs0': 11.100, 'logMs0Err': 0.018,
                              'beta': 0.470, 'betaErr': 0.008,
                              'delta': 0.393, 'deltaErr': 0.088,
                              'gamma': 2.510, 'gammaErr': 0.250,
                              'sigMs': 0.227, 'sigMsErr': 0.020,
                              'bCut': 2.46, 'bCutErr': 0.53,
                              'bSat': 8.72, 'bSatErr': 0.53,
                              'betaCut': 0.570, 'betaCutErr': 0.200,
                              'betaSat': 0.863, 'betaSatErr': 0.053}

            else:
                raise KeyError("Wrong SIG_MOD choice!!")

        else:
            raise KeyError("Wrong Redshift choice!!")

    def toMhalo(self, Ms):
        """Estimate halo mass via stellar mass."""
        param = self.param

        mRatio = Ms / (10.0 ** param['logMs0'])

        termB = np.log10(mRatio) * param['beta']
        termC = mRatio ** param['delta']
        termD = mRatio ** (param['gamma'] * -1.0) + 1.0

        logMh = param['logM1'] + termB + (termC / termD) - 0.50

        return logMh
