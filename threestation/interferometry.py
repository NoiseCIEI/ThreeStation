"""Perfrom interferometry (correlation or convolution)."""
import logging.config

import numpy as np
import scipy as sp
from scipy.fft import fft, ifft, fftfreq, next_fast_len

import pymodule as my
from pymodule.seis import geod_inv


logger = logging.getLogger(__name__)


def flip_nlag(tr1, tr2, lags):
    """
    Flip negative lags for correlation.
    """
    if lags == 'NN':
        tr1.data = tr1.data[::-1]
        tr2.data = tr2.data[::-1]
    elif lags == 'PN':
        tr2.data = tr2.data[::-1]
    elif lags == 'NP':
        tr1.data = tr1.data[::-1]

    return


def _cut(tr, lag, b, e):
    if lag == 'N':
        tmp = b
        b = - e
        e = - tmp

    return my.seis.sliced(tr, b=b, e=e)


def overlap(tr1, tr2, lags):
    _b1 = tr1.stats.sac.b
    _e1 = tr1.stats.sac.e
    _b2 = tr2.stats.sac.b
    _e2 = tr2.stats.sac.e

    b1 = min(abs(_b1), abs(_e1))
    b2 = min(abs(_b2), abs(_e2))
    # Use ceil & floor to ensure same length
    b = np.ceil(max(b1, b2))

    e1 = max(abs(_b1), abs(_e1))
    e2 = max(abs(_b2), abs(_e2))
    e = np.floor(min(e1, e2))

    if b1 != b2:
        tr1 = _cut(tr1, lags[0], b, e)
        tr2 = _cut(tr2, lags[1], b, e)

    return tr1, tr2


def _src_distaz(lon_r1, lat_r1, lon_r2, lat_r2, lon_s, lat_s):
    """
    Return source location relative to the center of a pair of receivers.
    """
    lon_c = (lon_r1 + lon_r2) / 2
    lat_c = (lat_r1 + lat_r2) / 2
    dist_cr, _, baz_cr = geod_inv(lat2=lat_c, lon2=lon_c, lat1=lat_r1, lon1=lon_r1)
    dist_cs, _, baz_cs = geod_inv(lat2=lat_c, lon2=lon_c, lat1=lat_s, lon1=lon_s)

    baz = (baz_cs - baz_cr) % 360

    return dist_cr, dist_cs, baz


def triangle_edges(lon_r1, lat_r1, lon_r2, lat_r2, lon_s, lat_s):
    """
    Return lengths of 3 sides of a triangle on sphere [km].
    """
    r1s, _, _ = geod_inv(lat1=lat_r1, lon1=lon_r1, lat2=lat_s, lon2=lon_s)
    r2s, _, _ = geod_inv(lat1=lat_r2, lon1=lon_r2, lat2=lat_s, lon2=lon_s)
    r12, _, _ = geod_inv(lat1=lat_r1, lon1=lon_r1, lat2=lat_r2, lon2=lon_r2)

    return r1s, r2s, r12


def stationary_phase_zone(lon_r1, lat_r1, lon_r2, lat_r2, lon_s, lat_s, **kwargs):
    """
    If a source-station is in stationary phase zone (SPZ).
    """
    op = kwargs.get('operator', 'correlation').lower()
    method = kwargs.get('method', 'distance').lower()
    max_drpct = kwargs.get('max_drpct', 1)
    max_dr = kwargs.get('max_dr', 30)
    max_deg = kwargs.get('max_deg', 8)
    return_srcdir = kwargs.get('return_srcdir', True)
    min_srdist = kwargs.get('min_srdist', 30)

    if op in ['corr', 'correlation']:
        use_corr = True
        use_conv = False
    elif op in ['conv', 'convolution']:
        use_corr = False
        use_conv = True
    else:
        raise ValueError(f'Unknown operator: {op}')

    r1s, r2s, r12 = triangle_edges(lon_r1, lat_r1, lon_r2, lat_r2, lon_s, lat_s)
    srdist = min(r1s, r2s)
    if srdist < min_srdist:
        return [False] + [None]*5
    dist_cr, dist_cs, _theta = \
        _src_distaz(lon_r1, lat_r1, lon_r2, lat_r2, lon_s, lat_s)
    op = op.lower()
    if use_corr:
        dr = abs(r1s - r2s) - r12
        theta = _theta
    elif use_conv:
        dr = (r1s + r2s) - r12
        theta = np.rad2deg(np.arccos(r12 / (r1s+r2s)))

    if method in ['dist', 'distance']:
        inspz = abs(dr) < max_drpct/100 * r12
    elif method in ['const', 'constant']:
        inspz = abs(dr) < max_dr
    elif method in ['az', 'azimuth']:
        if use_corr:
            inspz = False
            for k in [0, 180, 360]:
                if abs(theta-k) < max_deg:
                    inspz = True
                    break
        elif use_conv:
            inspz = theta < max_deg
    else:
        raise ValueError(f'Unknown SPZ method: {method}')

    if return_srcdir:
        if r1s < r2s:
            dir_src = 1
        else:
            dir_src = 2
        return inspz, dr, theta, r12, srdist, dir_src
    else:
        return inspz, dr, theta, r12, srdist


def _phase(freq, dr, per, pv, **kwargs):
    """
    Return phase (:math:`\delta\phi`) for phase shift.
    Be careful about omega and f.
    """
    kind = kwargs.get('kind', 'linear')
    fill_value = kwargs.get('fill_value', 'extrapolate')
    assume_sorted = kwargs.get('assume_sorted', True)

    f0 = 1 / per
    f0 = f0[::-1]
    pv = pv[::-1]
    fitp = sp.interpolate.interp1d(
        f0,
        pv,
        kind=kind,
        fill_value=fill_value,
        assume_sorted=assume_sorted,
    )
    c = fitp(freq)
    ph = 2*np.pi*freq * dr / c
    nfreq = freq.size
    if nfreq % 2 == 0:
        i = int(nfreq / 2)
        ph[i+1:] = -ph[1:i][::-1]
        ph[i] = 0
    else:
        i = int((nfreq+1) / 2)
        ph[i:] = -ph[1:i][::-1]

    return ph


def phase_shift(delta, dr, per, pv, spc=None, data=None):
    """
    Perfrom phase shift.

    :param delta: Sample spacing.
    :param dr: Difference of differential or sum of source-receiver distances
        and receiver-distance.
    """
    if spc is None:
        nfreq = next_fast_len(data.size)
        spc = fft(data, nfreq)
    else:
        nfreq = spc.size
    freq = fftfreq(nfreq, d=delta)
    dph = _phase(freq, dr, per, pv)
    spc = spc * np.exp(1j * dph)
    data_ps = ifft(spc, nfreq).real

    return data_ps


def xc_ps(tr1, tr2, delta, dr, per, pv, **kwargs):
    """
    Perform correlation of two :class:`~obspy.core.trace.Trace`,
    then apply phase shift in frequency domain and transform back.

    https://docs.obspy.org/packages/autogen/obspy.signal.cross_correlation.correlate.html#obspy.signal.cross_correlation.correlate
    """
    demean = kwargs.get('demean', False)
    op = kwargs.get('operator', 'corr').lower()

    delta = tr1.stats.delta
    npts = min(tr1.stats.npts, tr2.stats.npts)
    shift = kwargs.get('shift', npts - 1)

    mid = (tr1.stats.npts + tr2.stats.npts - 1) // 2
    if shift > mid:
        raise ValueError('Such a large shift is not possible without zero padding')

    # Interchange to make convention same as SAC
    a1 = tr2.data
    if op in ['corr', 'correlation']:
        a2 = tr1.data
    elif op in ['conv', 'convolution']:
        a2 = tr1.data[::-1]

    if demean:
        a1 = a1 - np.mean(a1)
        a2 = a2 - np.mean(a2)

    spc = my.signal.conv_spc(a1, a2[::-1])
    xc = phase_shift(delta, dr, per, pv, spc=spc)

    return xc[mid-shift : mid+shift + xc.size%2]


def pick_lag(a, dir_src):
    """
    Pick positive or negative lag of correlation based on source direction.
    """
    mid = int((a.size-1) / 2)
    if dir_src == 1:
        return a[mid:]
    else:
        return a[:mid+1][::-1]
