"""Preprocess two-station interferograms (:math:`I_2`)."""
import logging.config
# from sys import exit

import bottleneck as bn
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import obspy
import scipy as sp
from scipy.fft import fft, ifft, fftfreq, next_fast_len


logger = logging.getLogger(__name__)


def one_bit(a):
    """
    Return sign values of an array, i.e., one-bit normalization.
    """
    return np.sign(a)


def whiten(tr, Tmin=1, Tmax=150, freq_width=.0004, brute=True,
           frac=.2, epsilon=1e-8, **kwargs):
    """
    Spectral whitening using running absolute mean in frequency domain.

    https://github.com/NoiseCIEI/Seed2Cor/blob/master/src/Whiten.c
    https://github.com/bgoutorbe/seismic-noise-tomography/blob/8f26ff827bee8a411038e33d93b59bffbc88c0a7/pysismo/pscrosscorr.py#L306
    https://www.mathworks.com/matlabcentral/fileexchange/65345-spectral-whitening
    https://github.com/ROBelgium/MSNoise/blob/27749e2914b30b1ab8278054645677444f10e309/msnoise/move2obspy.py#L135

    :param freq_width: length of averaging window in frequency domain
    :param epsilon: minimum value to avoid zero division
    :param return_spc: if return spectra for plot
    """
    npts = tr.stats.npts
    sr = tr.stats.sampling_rate
    dom = sr / npts
    winlen = int(round(freq_width / dom))

    nfreq = next_fast_len(npts)
    spc = fft(tr.data, nfreq)
    spc_am = np.abs(spc)
    spc_ph = np.unwrap(np.angle(spc))

    if brute:
        weight = spc_am
        spc_w = np.exp(1j * spc_ph)
    elif winlen < 2:
        weight = spc_am
        spc_w = np.exp(1j * spc_ph)
        logger.debug('Window too short')
    else:
        weight = bn.move_mean(spc_am, winlen, 1)
        ind = weight < epsilon
        weight[ind] = epsilon
        spc_w = spc / weight
        spc_w[ind] = 0

    f2 = 1 / Tmax
    f3 = 1 / Tmin
    f1 = f2 * (1-frac)
    f4 = f3 * (1+frac)
    if f4 > sr / 2:
        logger.warning('Whiten band upper end out of range! Corrected to Nyquist.')
        f4 = sr / 2
        if f3 >= f4:
            f3 = f4 * (1-frac)
            if f3 < f2:
                f3 = f2
    freqs = fftfreq(nfreq, d=tr.stats.delta)
    spc_w *= obspy.signal.invsim.cosine_taper(
            npts=nfreq,
            freqs=freqs,
            flimit=[f1, f2, f3, f4],
        )

    # Hermitian symmetry (because the input is real)
    spc_w[-(nfreq//2)+1:] = spc_w[1:(nfreq//2)].conjugate()[::-1]

    whitened = ifft(spc_w, nfreq)[:npts]

    if kwargs.get('plot', False):
        xlim = [0, .2]
        npts = tr.stats.npts
        delta = tr.stats.delta
        t = np.arange(0, npts * delta, delta)

        fig = plt.figure(figsize=(12, 8))
        gs = mpl.gridspec.GridSpec(3, 1)
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[2, 0], sharex=ax2)

        tr.filter('bandpass', freqmin=f2, freqmax=f3, zerophase=True)
        ax1.plot(t, tr.normalize().data, c='k', label="Raw", ls='--', lw=1)
        ax1.plot(t, whitened/np.abs(whitened).max(), c='r', label="Whitened", lw=1)
        ax1.set_xlabel("Time (s)")
        ax1.set_xlim(0, tr.stats.sac.dist)
        ax1.set_ylim(-1, 1)

        ifreq = np.argsort(freqs)
        ax2.plot(freqs[ifreq], spc_am[ifreq],
                 alpha=.5, lw=1, c='gray', ls='--', label="Raw")
        ax2.plot(freqs[ifreq], weight[ifreq],
                 alpha=.5, lw=2, c='g', label="Weight")
        ax2.set_xlim(xlim[0], xlim[1])

        ax3.plot(freqs[ifreq], np.abs(spc_w)[ifreq], lw=1, c='b', label="Whitened")
        for ax in [ax2, ax3]:
            ax.set_ylim(0)
            ax.set_xlabel("Frequency (Hz)")
            for per in [8, 16, 26]:
                ax.axvline(1/per, ls='--', alpha=.5, c='k')
        for ax in [ax1, ax2, ax3]:
            ax.legend(loc='upper right')

        pair = f'{tr.stats.sac.kevnm.strip()}_{tr.stats.sac.kstnm.strip()}'
        fig.suptitle(pair, y=1.02)
        # fig.savefig(f'/work2/szhang/US/Plot/fig4exp/US/sw_{pair}.pdf')
        plt.show()

    tr.data = whitened
    tr.taper(max_percentage=0.05, type='hann', side='both')

    return tr


def mute(a, t1, t2, sr, precursor_only=True):
    """
    Mute outiside [t1, t2].

    :param sr: sampling rate
    :param precursor_only: Only mute precursory (< t1) for caculating SNR.
    """
    n1 = int(t1 * sr)
    n2 = int(t2 * sr)
    a[:n1] = 0
    if not precursor_only:
        a[n2:] = 0

    return a
