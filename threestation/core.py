"""Core functions for three-station interferometry."""
from os.path import basename, exists, dirname, join
import multiprocessing as mp
import functools
import itertools as it
import logging.config
from sys import exit

import numpy as np
from obspy.core import (
    read,
    Stream,
    Trace,
)
from tqdm import tqdm

import pymodule as my
from .config import (
    # Constants
    CONV,
    CORR,
    DEF_SHD,
    DIROUT,
    KEY2SHD,
    PARAM,
    META,
    RECEIVER_STATION,
    SOURCE_STATION,
    USE_CW,
    USE_DW,

    # Functions
    get_fnm,
    get_pred_pv,
)
from .preprocess import (
    one_bit,
    mute,
    whiten,
)
from .interferometry import (
    flip_nlag,
    stationary_phase_zone,
    triangle_edges,
    xc_ps,
    phase_shift,
    pick_lag,
    overlap,
)


logger = logging.getLogger(__name__)


# To be filled
STNM2SRC = {}
STNM2I2 = {}
DEST2LAG = {}
LST_LAG = list(DEST2LAG.keys())
PAIR2SRC = {}


def _receiver_station_pairs():
    """
    Return all possible combinations of two receiver-stations.
    """
    if 'group' in RECEIVER_STATION.columns:
        return it.product(*[list(i['net_sta']) for _, i in RECEIVER_STATION.groupby('group')])
    else:
        return it.combinations(list(RECEIVER_STATION['net_sta']), 2)


def get_two_station_interferogram():
    """
    Get two-station ambient noise interferograms, :math:`I^{AN}_2`.
    `STNM2SRC` will map station names to a set of stations
    with which the station have correlation.
    `STNM2I2` will map station names to paths to I2 files.
    `STNM2I2` can be saved to the metadata directory.
    """
    global STNM2SRC
    global STNM2I2

    if PARAM['skip']['find_I2']:
        logger.debug('Skip getting I2')
        return

    logger.info('# Get two-station interferograms for receiver-stations.')
    STNM2SRC.update({rec: set() for rec in RECEIVER_STATION['net_sta']})
    STNM2I2.update({rec: set() for rec in RECEIVER_STATION['net_sta']})
    for receiver, source in tqdm(list(it.product(list(RECEIVER_STATION['net_sta']), list(SOURCE_STATION['net_sta'])))):
        # Avoid using auto-correlations
        if receiver == source:
            continue
        raw = get_fnm('I2', receiver, source)
        if exists(raw):
            STNM2I2[receiver].add(raw)
            STNM2SRC[receiver].add(source)

    if PARAM['write']['meta']:
        my.fio.wpk(get_fnm('I2_path'), STNM2I2)

    return


def _cut_ends(I2):
    """
    Return ends to cut.
    """
    p = PARAM['cut']

    dist = I2.stats.sac.dist
    if USE_DW:
        bp = max(0, dist/p['vmax'] - p['bfact_dw']*p['Tmax'])
        ep = dist/p['vmin'] + p['efact_dw']*p['Tmax']
    else:
        bp = dist/p['vmin'] + p['efact_cw']*p['Tmax']
        ep = p['te']

    return bp, ep, - bp, - ep


def _cut_signal(trnm, dests, epsilon=1e-4):
    """
    Return cut of positive & negative lags.
    """
    I2 = read(trnm, format='SAC')[0]

    # Check window length
    sr = I2.stats.sampling_rate
    bp, ep, en, bn = _cut_ends(I2)
    tlen = ep - bp
    if USE_CW and (tlen < PARAM['cut']['min_len']):
        logger.debug(f'{basename(trnm)}: tlen = {tlen:.1f} s')
        return {}

    # SNR threshold on I2
    sym = my.seis.sym_xc(I2)
    try:
        _snr = I2.stats.sac[KEY2SHD['snr']]
    except (KeyError, AttributeError):
        _snr = my.seis.snr(sym, **PARAM['snr'])
        I2.stats.sac[KEY2SHD['snr']] = _snr
        I2.write(trnm)

    if _snr < PARAM['snr']['min']:
        logger.debug(f'{basename(trnm)}: SNR = {_snr:.1f}')
        return {}

    # Check delta
    _delta = I2.stats.delta
    delta2 = PARAM['cut']['delta']
    if abs(_delta - delta2) > epsilon:
        logger.warning(f'Resample {trnm} (delta={_delta:g})')
        I2.resample(sampling_rate=1/delta2)

    # Cut
    if USE_CW:
        sym.data = sym.data[int(bp*sr):]
        sym.stats.sac.b = bp
        sym.taper(max_percentage=0.05, type='hann', side='both')

    if len(dests) == 1:
        lags = [sym]
    elif len(dests) == 2:
        e = I2.stats.sac.e
        plag = my.seis.sliced(I2, 0, e)
        plag.stats.sac.b = 0
        plag.stats.sac.e = e
        nlag = my.seis.sliced(I2, -e, 0)
        nlag.stats.sac.b = -e
        nlag.stats.sac.e = 0
        lags = [plag, nlag]

    # 1-bit
    if PARAM['preproc'].get('onebit', False):
        for lag in lags:
            lag.data = one_bit(lag.data)

    # Whiten
    if PARAM['preproc']['whiten'].get('val', False):
        for lag in lags:
            lag = whiten(lag, **PARAM['preproc']['whiten'])

    # Mute
    if USE_DW and getattr(PARAM['cut'], 'mute', True):
        kwargs_mute = {
            't1': bp,
            't2': ep,
            'sr': sr,
            'precursor_only': PARAM['cut']['mute_prc_only'],
        }
        for i, lag in enumerate(lags):
            # Negative lag
            if i == 1:
                lag.data = mute(lag.data[::-1], **kwargs_mute)[::-1]
            else:
                lag.data = mute(lag.data, **kwargs_mute)

    d2l = {}
    for lag, dest in zip(lags, dests):
        d2l[dest] = lag

    return d2l


def cut_signal():
    """
    Cut signal from :math:`I_2`.
    `DEST2LAG` is a dictionary mapping filenames of :math:`I_2` to corresponding traces.
    This dictionary saves all :math:`I_2` traces for later use, to avoid reading
    the same trace twice and to ease multiprocessing.
    The chosen signals can be saved to the path specified in
    :func:`threestation.config.fnm`.

    """
    global DEST2LAG

    if PARAM['skip']['cut']:
        logger.debug('Skip cutting signal')
        return

    logger.info('# Cut signal from two-station interferograms')

    args_mp = []
    for receiver in RECEIVER_STATION['net_sta']:
        if PARAM['write']['lag']:
            my.sys_tool.mkdir(join(DIROUT, receiver))
        for I2 in STNM2I2[receiver]:
            dests = get_fnm('I2_lag_raw', receiver, I2=I2)
            # Avoid duplication when sta2 is also a receiver
            sta2 = basename(dirname(I2))
            if (
                (sta2 != receiver)
                and (sta2 in RECEIVER_STATION['net_sta'])
                and (I2 in STNM2I2[sta2])
            ):
                continue
            args_mp.append([I2, dests])

    with mp.Pool(PARAM['misc']['ncpu']) as p:
        _dest = p.starmap(_cut_signal, tqdm(args_mp))

    for i in _dest:
        if i:
            DEST2LAG.update(i)

    if PARAM['write']['lag']:
        for dest, lag in DEST2LAG.items():
            lag.write(dest)

    return


def find_common_source():
    """
    Find common source-stations for receiver-stations.
    `PAIR2SRC` will map names of station-pair to the corresponding
    set of source-stations, which can be saved in the metadata folder.
    """
    global PAIR2SRC

    if PARAM['skip']['find_source']:
        logger.debug('Skip finding common source-stations')
        return

    logger.info('# Find common source-stations for receiver-station pairs')
    for sta1, sta2 in tqdm(list(_receiver_station_pairs())):
        PAIR2SRC[f'{sta1}_{sta2}'] = STNM2SRC[sta1] & STNM2SRC[sta2]

    if PARAM['write']['meta']:
        my.fio.wpk(get_fnm('source-station'), PAIR2SRC)

    return


def _source_specifc_interferogram(trnm1, trnm2, rec1, rec2, src, lags, **kwargs):
    """
    Construct a source-specific interferogram (:math:`C_3`), i.e.,
    interferogram of two :math:`I_2`.

    :param rec1: Name of the first receiver-station.
    :param rec2: Name of the second receiver-station.
    :param src: Name of the source-station.
    :param lags: Use which lags for I3.
    :param dir_src: If save source direction to SAC header.
    """
    dest = get_fnm('C3', rec1, sta2=rec2, sta3=src, lags=lags)

    if PARAM['skip']['C3'] and exists(dest):
        logger.debug(f'{dest} already exists.')
        if PARAM['write']['stack']:
            return dest, read(dest, format='SAC')[0]
        else:
            return None, None

    tr1 = DEST2LAG[trnm1]
    tr2 = DEST2LAG[trnm2]

    # Find common part
    if USE_CW:
        tr1, tr2 = overlap(tr1, tr2, lags)

    # Flip negative lag
    if PARAM['interferometry']['flip_nlag']:
        flip_nlag(tr1, tr2, lags)

    # Do interferometry
    if USE_DW and PARAM['interferometry']['phase_shift']:
        kwa_ps = {
            'delta': tr1.stats.delta,
            'dr': kwargs.get('dr'),
            'per': kwargs.get('phprper'),
            'pv': kwargs.get('phprvel'),
        }
        if CONV:
            C3 = xc_ps(tr1, tr2, **kwa_ps, **PARAM['interferometry'])
        elif CORR:
            xc = my.seis.x_cor(tr1, tr2, **PARAM['interferometry'])
            C3 = pick_lag(xc, kwargs.get('dir_src'))
            C3 = phase_shift(data=C3, **kwa_ps)
    else:
        C3 = my.seis.x_cor(tr1, tr2, **PARAM['interferometry'])
        if USE_DW and CORR and PARAM['interferometry']['pick_lag']:
            C3 = pick_lag(C3, kwargs.get('dir_src'))

    # Make header
    b = - int(np.floor(C3.size/2) / tr1.stats.delta)
    if USE_CW:
        if PARAM['interferometry'].get('Welch', False):
            b = - PARAM['interferometry']['subwin']
    if USE_DW:
        if CONV or PARAM['interferometry']['pick_lag']:
            b = 0
    if PARAM['interferometry']['symmetric'] or CONV:
        nsided = 1
    else:
        nsided = 2

    header = my.seis.sachd(**{
        'b': b,
        'e': int(b + C3.size*tr1.stats.delta),
        'delta': tr1.stats.delta,
        'npts': C3.size,
        'kevnm': rec1,
        'evlo': META[META['net_sta'] == rec1]['lon'].iloc[0],
        'evla': META[META['net_sta'] == rec1]['lat'].iloc[0],
        'knetwk': rec2.split('_')[0],
        'kstnm': rec2,
        'stlo': META[META['net_sta'] == rec2]['lon'].iloc[0],
        'stla': META[META['net_sta'] == rec2]['lat'].iloc[0],
        'dist': kwargs.get('dist', DEF_SHD),
        KEY2SHD['net_rec']: rec2.split('_')[0],  # to be consistent with I2
        KEY2SHD['src_sta']: src.split('_')[1],
        KEY2SHD['src_net']: src.split('_')[0],
        KEY2SHD['nsided']: nsided,
        KEY2SHD['dr']: kwargs.get('dr', DEF_SHD),
        KEY2SHD['theta']: kwargs.get('theta', DEF_SHD),
        KEY2SHD['dir_src']: kwargs.get('dir_src', DEF_SHD),
        KEY2SHD['min_srdist']: kwargs.get('min_srdist', DEF_SHD),
        })

    # Make Trace
    C3_tr = Trace(header=header, data=C3)
    if USE_CW and PARAM['interferometry']['symmetric']:
        C3_tr = my.seis.sym_xc(C3_tr)

    if CONV and PARAM['interferometry'].get('trim_conv', True):
        i = int(PARAM['cut']['te'] / PARAM['cut']['delta']) + 1
        C3_tr.data = C3_tr.data[:i]

    return dest, C3_tr


def _source_specifc_interferogram_pair(rec1, rec2, src, phprper=None, phprvel=None):
    """
    Compute source-specific interferograms (:math:`C_3`) for a receiver-pair
    using 2 or 4 lags of two-station interferograms (:math:`I_2`).

    :param rec1: Name of the first receiver-station.
    :param rec2: Name of the second receiver-station.
    :param src: Name of the source-station.
    :param phprper: Periods of prior dispersion for phase shift.
    :param phprvel: Phase velocities of prior dispersion for phase shift.
    """
    dist = None
    dr = DEF_SHD
    theta = DEF_SHD
    dir_src = DEF_SHD
    min_srdist = DEF_SHD

    lon_1 = META[META['net_sta'] == rec1]['lon'].iloc[0]
    lat_1 = META[META['net_sta'] == rec1]['lat'].iloc[0]
    lon_2 = META[META['net_sta'] == rec2]['lon'].iloc[0]
    lat_2 = META[META['net_sta'] == rec2]['lat'].iloc[0]
    lon_s = META[META['net_sta'] == src]['lon'].iloc[0]
    lat_s = META[META['net_sta'] == src]['lat'].iloc[0]
    res = stationary_phase_zone(
        lon_1, lat_1, lon_2, lat_2, lon_s, lat_s,
        **PARAM['interferometry']
    )
    inspz, dr, theta, dist, min_srdist = res[:5]
    if PARAM['interferometry'].get('dir_src', True):
        dir_src = res[5]

    # Limit source station in SPZ
    if PARAM['interferometry'].get('spz', False):
        if not inspz:
            return [], []

    if USE_CW and getattr(PARAM['cut'], 'restrict_src_dist', True):
        r1s, r2s, dist = triangle_edges(
            lon_1, lat_1, lon_2, lat_2, lon_s, lat_s,
        )
        for _r in [r1s, r2s]:
            t_coda = (_r/PARAM['cut']['vmin']
                      + PARAM['cut']['efact_cw'] * PARAM['cut']['Tmax'])
            t_snr = dist/PARAM['snr']['vmin'] + PARAM['cut']['cw4snr']
            if t_snr > (PARAM['cut']['te'] - t_coda):
                logger.debug(f'{src} too far from {rec1}-{rec2}')
                return [], []

    if dist is None:
        dist, _, _ = my.seis.geod_inv(lon_1, lat_1, lon_2, lat_2)

    interferogram = functools.partial(
        _source_specifc_interferogram,
        rec1=rec1,
        rec2=rec2,
        src=src,
        dr=dr,
        theta=theta,
        dist=dist,
        dir_src=dir_src,
        min_srdist=min_srdist,
        phprper=phprper,
        phprvel=phprvel,
    )

    dest = [None]
    I3 = [None]

    if PARAM['interferometry']['nlag'] == 4:
        s1 = get_fnm('I2_lag_proc', rec1, sta2=src, pre=PARAM['pfx']['sym'])
        s2 = get_fnm('I2_lag_proc', rec2, sta2=src, pre=PARAM['pfx']['sym'])
        if (s1 in DEST2LAG) and (s2 in DEST2LAG):
            dest, ss = interferogram(trnm1=s1, trnm2=s2, lags='SS')
            dest = [dest]
            I3 = [ss]

    elif PARAM['interferometry']['nlag'] == 2:
        p1 = get_fnm('I2_lag_proc', rec1, sta2=src, pre=PARAM['pfx']['plag'])
        p2 = get_fnm('I2_lag_proc', rec2, sta2=src, pre=PARAM['pfx']['plag'])
        n1 = get_fnm('I2_lag_proc', rec1, sta2=src, pre=PARAM['pfx']['nlag'])
        n2 = get_fnm('I2_lag_proc', rec2, sta2=src, pre=PARAM['pfx']['nlag'])
        if (
            (p1 in DEST2LAG)
            and (p2 in DEST2LAG)
            and (n1 in DEST2LAG)
            and (n2 in DEST2LAG)
        ):
            dest_p, pp = interferogram(trnm1=p1, trnm2=p2, lags='PP')
            dest_n, nn = interferogram(trnm1=n1, trnm2=n2, lags='NN')
            dest = [dest_p, dest_n]
            I3 = [pp, nn]

    else:
        raise ValueError(f"Unkown nlag {PARAM['interferometry']['nlag']} (2 or 4)")

    if None not in dest:
        return dest, I3
    else:
        return [], []


def _three_station_interferometry_pair(rec1, rec2):
    """
    Calculate source-specific interferogram (:math:`C_3`)
    and stack for :math:`I_3` for a station-pair.
    """
    dest = get_fnm('I3', rec1, sta2=rec2)
    if PARAM['skip']['I3'] and exists(dest):
        logger.debug(f'Skip {dest}')
        return dest

    phprper = None
    phprvel = None
    if USE_DW and PARAM['interferometry']['phase_shift']:
        phprper, phprvel = get_pred_pv(rec1, rec2)

    dir_src = PARAM['interferometry'].get('dir_src', True)
    dest_I3s = []
    I3s = Stream()
    nsrc = 0
    n1 = 0
    n2 = 0
    for src in PAIR2SRC[f'{rec1}_{rec2}']:
        dest_I3, I3 = _source_specifc_interferogram_pair(
            rec1, rec2, src,
            phprper=phprper, phprvel=phprvel,
        )
        if I3:
            dest_I3s.extend(dest_I3)
            I3s.extend(I3)
            nsrc += 1
            if dir_src:
                if I3[0].stats.sac[KEY2SHD['dir_src']] == 1:
                    n1 += 1
                else:
                    n2 += 1

    # Do stack
    if PARAM['write']['stack']:
        my.sys_tool.mkdir(join(DIROUT, PARAM['dir']['I3'], rec1))
        if nsrc < PARAM['stack']['min_src']:
            logger.debug(f'Insufficient ({nsrc}) source-stations for {rec1}-{rec2}')
            return

        logger.debug(f'{nsrc} source-stations for {rec1}-{rec2}')

        stats = I3s[0].stats.copy()
        stats.sac[KEY2SHD['nsrc']] = nsrc
        if dir_src:
            stats.sac[KEY2SHD['nsrc_dir1']] = n1
            stats.sac[KEY2SHD['nsrc_dir2']] = n2
        for key in [
            'src_net',
            'src_sta',
            'dir_src',
            'dr',
            'theta',
        ]:
            stats.sac.pop(KEY2SHD[key])

        sk = my.seis.stack(
            I3s,
            stats=stats,
            kw_stack={**PARAM['stack'], **KEY2SHD},
            kw_snr=PARAM['snr'],
        )
        if sk is None:
            logger.warning(f'{rec1}-{rec2}: Failed for SNR')
            return
        else:
            sk.write(dest, format='SAC')

        if PARAM['stack']['rand']:
            my.sys_tool.mkdir(join(
                DIROUT,
                PARAM['dir']['I3_rand'],
                rec1,
                rec2,
            ))
            dest_rand = get_fnm('I3_rand', rec1, sta2=rec2)
            for i, _sk in enumerate(my.seis.rand_stack(
                I3s,
                stats=stats,
                kw_stack={**PARAM['stack'], **KEY2SHD},
                kw_snr=PARAM['snr'],
            )):
                _sk.write(f'{dest_rand}.{i}', format='SAC')

    # To save SNR in header after stack
    if PARAM['write']['C3']:
        my.sys_tool.mkdir(join(DIROUT, rec1, rec2))
        for nm, tr in zip(dest_I3s, I3s):
            tr.write(nm, format='SAC')

    return dest


def three_station_interferometry():
    """
    Calculate source-specific interferogram (:math:`C_3`)
    and stack for :math:`I_3` for all station pairs.
    """
    logger.info('# Construct and stack source-specific three-station interferograms')
    with mp.Pool(PARAM['misc']['ncpu']) as p:
        p.starmap(
            _three_station_interferometry_pair,
            tqdm(list(_receiver_station_pairs())),
        )

    return


def main():
    """
    Workflow of three-station interferometry:

    1. Find paths to two-station interferograms (:math:`I_2`).
    #. Cut and preprocess signals from :math:`I_2`.
    #. Find common source-stations for each pair of receiver-stations.
    #. Calculate source-specific interferogram (:math:`C_3`)
        and stack for :math:`I_3`.

    """
    get_two_station_interferogram()
    cut_signal()
    find_common_source()
    three_station_interferometry()

    return
