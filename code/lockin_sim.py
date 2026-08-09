"""
lockin_sim.py
-------------
Is a lock-in amplifier worth buying for our high-pressure CW-ODMR setup?

The honest version of that question is: what noise currently limits us?  A
lock-in cannot beat shot noise -- it can only move the measurement band above
the technical 1/f noise (laser intensity drift, focus/DAC thermal drift,
background wander).  So the decision is set by one number, the knee frequency
f_knee of the photodetector noise spectrum, compared with how fast we can
modulate.

This module answers it by simulation, with no experimental input:

  * technical (multiplicative) noise with one-sided relative PSD
        S(f) = S0 (1 + f_knee/f),
    S0 pinned to the shot-noise floor so that the setup is shot-noise limited
    at high frequency and f_knee is the ONLY free knob;
  * shot noise added independently per sample;
  * three detection schemes, all given the SAME total measurement time and the
    SAME noise realisation, so the comparison is fair:
      DC     : slow frequency sweep, MW always on, Lorentzian fit
      ON-OFF : MW chopped at f_mod, square-wave demodulation, Lorentzian fit
      FM     : MW frequency dithered at f_mod, sine demodulation, dispersive
               (derivative) fit
  * figure of merit = the standard deviation of the fitted line centre D over
    many noise realisations, i.e. exactly the quantity that limits a pressure
    or magnetic-field reading.

Everything is phenomenological only in the noise model; the estimators are
the ones an experiment would actually use.
"""

import numpy as np
from scipy.optimize import curve_fit

# ------------------------------------------------------------- defaults ----
FS = 20_000.0          # sample rate of the photodetector record (Hz)
T_TOTAL = 0.5          # total measurement time per ODMR sweep (s)
N_PTS = 51             # frequency points in the sweep
R_PHOT = 1.0e8         # detected photon rate (counts/s)
CONTRAST = 0.03        # ODMR contrast (3% -- realistic at 120 GPa)
LINEWIDTH = 8.0        # ODMR FWHM (MHz)
D_TRUE = 2870.0        # line centre (MHz)
SPAN = 3.0             # sweep half-width in units of the linewidth
F_MOD = 2000.0         # modulation frequency (Hz)
FM_DEPTH = 0.5         # FM dither amplitude in units of the linewidth (near-optimal)


# ----------------------------------------------------------- noise model ---
def technical_noise(n_samp, fs, S0, f_knee, rng):
    """Relative intensity noise with one-sided PSD S0 (1 + f_knee/f).

    Built by shaping white noise in the frequency domain.  The DC bin is set
    to zero: a constant offset is absorbed by the baseline term of every fit
    and is not what distinguishes the schemes.
    """
    f = np.fft.rfftfreq(n_samp, 1.0 / fs)
    S = np.zeros_like(f)
    S[1:] = S0 * (1.0 + f_knee / f[1:])
    amp = np.sqrt(S * fs * n_samp / 4.0)
    X = amp * (rng.normal(size=f.size) + 1j * rng.normal(size=f.size))
    X[0] = 0.0
    return np.fft.irfft(X, n_samp)


def shot_floor_psd(R):
    """One-sided relative PSD of shot noise on a photon rate R: 2/R per Hz."""
    return 2.0 / R


# ------------------------------------------------------------- lineshape ---
def lorentz(nu, nu0, w):
    return 1.0 / (1.0 + ((nu - nu0) / (0.5 * w)) ** 2)


def _fit_dip(nu, y, w, positive):
    """Fit Lorentzian (dip or peak) + linear baseline; return the centre.

    The linewidth is held fixed: in practice it is known from a separate
    high-SNR calibration scan, and freeing it only adds fit noise that is
    common to all three schemes.
    """
    mid = nu.mean()
    sgn = 1.0 if positive else -1.0

    def model(v, amp, c, b0, b1):
        return amp * lorentz(v, c, w) + b0 + b1 * (v - mid)

    p0 = [sgn * (y.max() - y.min()), nu[np.argmax(sgn * y)], np.median(y), 0.0]
    try:
        p, _ = curve_fit(model, nu, y, p0=p0, maxfev=2000)
    except Exception:
        return np.nan
    return p[1] if nu[0] < p[1] < nu[-1] else np.nan


_FM_PHASE = np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False)
_FM_SIN = np.sin(_FM_PHASE)


def fm_harmonic(nu, c, w, depth):
    """Exact first-harmonic response of a Lorentzian dip to a sine dither.

    Using the analytic derivative instead is only valid for depth << w; at the
    depths that maximise the signal the demodulated shape is a broadened
    pseudo-derivative, and fitting it with a pure derivative biases the centre.
    """
    x = np.atleast_1d(nu)[:, None] - c + depth * w * _FM_SIN[None, :]
    return -2.0 * np.mean(lorentz(x, 0.0, w) * _FM_SIN[None, :], axis=1)


def _fit_dispersive(nu, y, w, depth):
    """Fit the exact demodulated FM lineshape + linear baseline; return the centre."""
    mid = nu.mean()
    shape = fm_harmonic(nu, mid, w, depth)
    pp = shape.max() - shape.min()

    def model(v, amp, c, b0, b1):
        return amp * fm_harmonic(v, c, w, depth) + b0 + b1 * (v - mid)

    yd = y - np.median(y)
    c0 = 0.5 * (nu[np.argmax(yd)] + nu[np.argmin(yd)])
    p0 = [(y.max() - y.min()) / max(pp, 1e-12), c0, np.median(y), 0.0]
    try:
        p, _ = curve_fit(model, nu, y, p0=p0, maxfev=2000)
    except Exception:
        return np.nan
    return p[1] if nu[0] < p[1] < nu[-1] else np.nan


# ----------------------------------------------------------- timing grid ---
def dwell_samples(fs, T, n_pts, f_mod):
    """Samples per frequency point `m` and per modulation cycle `p`.

    Demodulation only rejects the DC level exactly if each dwell contains a
    WHOLE number of modulation cycles; otherwise the leftover partial cycle
    leaks the (large) mean PL into the demodulated output and the comparison
    between schemes becomes an artefact of the timing grid.  So `m` is snapped
    down to a multiple of `p`, and the modulation frequency to fs/p.
    """
    m0 = max(1, int(round(fs * T)) // n_pts)
    p = int(max(2, round(fs / f_mod)))
    p += p % 2                                # even -> exact 50% duty cycle
    m = max(p, (m0 // p) * p)
    return m, p


def cycles_per_dwell(fs=FS, T=T_TOTAL, n_pts=N_PTS, f_mod=F_MOD):
    m, p = dwell_samples(fs, T, n_pts, f_mod)
    return m // p


# --------------------------------------------------------- one experiment --
def simulate_sweep(f_knee, rng, scheme, tech_gain=1.0,
                   fs=FS, T=T_TOTAL, n_pts=N_PTS, R=R_PHOT, C=CONTRAST,
                   w=LINEWIDTH, D=D_TRUE, span=SPAN, f_mod=F_MOD,
                   fm_depth=FM_DEPTH, noise=None, shot=True):
    """Run one ODMR sweep with the given detection scheme; return (nu, signal).

    `noise` lets the caller reuse one technical-noise realisation across
    schemes so that the comparison is paired.
    """
    m, p = dwell_samples(fs, T, n_pts, f_mod)
    f_mod = fs / p                            # snapped to an integer sample count
    n_samp = m * n_pts
    t = np.arange(n_samp) / fs
    nu = D + np.linspace(-span * w, span * w, n_pts)

    if noise is None:
        S0 = tech_gain * shot_floor_psd(R)
        noise = technical_noise(n_samp, fs, S0, f_knee, rng)
    nu_t = np.repeat(nu, m)                   # swept frequency vs time

    # exact square wave / sine from the sample index: with an integer number of
    # cycles per dwell both references have exactly zero mean over each dwell,
    # so demodulation cancels the (large) mean PL identically.
    phase = np.arange(n_samp) % p
    square = np.where(phase < p // 2, 1.0, -1.0)
    sine = np.sin(2 * np.pi * phase / p)

    if scheme == 'dc':
        pl = R * (1.0 - C * lorentz(nu_t, D, w))
    elif scheme == 'onoff':
        pl = R * (1.0 - C * lorentz(nu_t, D, w) * (square > 0))
    elif scheme == 'fm':
        pl = R * (1.0 - C * lorentz(nu_t + fm_depth * w * sine, D, w))
    else:
        raise ValueError(scheme)

    y = pl * (1.0 + noise)
    if shot:
        y = y + rng.normal(scale=np.sqrt(np.maximum(pl, 1.0) * fs), size=n_samp)
    Y = y.reshape(n_pts, m)

    if scheme == 'dc':
        return nu, Y.mean(axis=1)
    if scheme == 'onoff':
        return nu, -(Y * square.reshape(n_pts, m)).mean(axis=1)   # positive peak
    return nu, (Y * sine.reshape(n_pts, m)).mean(axis=1)


def estimate_D(f_knee, rng, scheme, **kw):
    """One realisation -> one estimate of the line centre (MHz)."""
    nu, sig = simulate_sweep(f_knee, rng, scheme, **kw)
    w = kw.get('w', LINEWIDTH)
    if scheme == 'fm':
        return _fit_dispersive(nu, sig, w, kw.get('fm_depth', FM_DEPTH))
    return _fit_dip(nu, sig, w, positive=(scheme == 'onoff'))


def sigma_D(f_knee, n_mc=200, seed=0, schemes=('dc', 'onoff', 'fm'), **kw):
    """Monte-Carlo standard deviation of the fitted centre, per scheme.

    All schemes see the SAME technical-noise realisation in each trial.
    """
    rng = np.random.default_rng(seed)
    fs = kw.get('fs', FS)
    T = kw.get('T', T_TOTAL)
    n_pts = kw.get('n_pts', N_PTS)
    R = kw.get('R', R_PHOT)
    S0 = kw.get('tech_gain', 1.0) * shot_floor_psd(R)
    n_samp = dwell_samples(fs, T, n_pts, kw.get('f_mod', F_MOD))[0] * n_pts

    out = {s: [] for s in schemes}
    for _ in range(n_mc):
        noise = technical_noise(n_samp, fs, S0, f_knee, rng)
        for s in schemes:
            out[s].append(estimate_D(f_knee, rng, s, noise=noise, **kw))
    return {s: float(np.nanstd(np.array(v) - kw.get('D', D_TRUE)))
            for s, v in out.items()}


if __name__ == '__main__':
    print(f'{"f_knee[Hz]":>10} {"DC":>10} {"ON-OFF":>10} {"FM":>10}   (sigma_D in MHz)')
    for fk in (0.0, 1.0, 10.0, 100.0, 1000.0):
        s = sigma_D(fk, n_mc=120, seed=1)
        print(f'{fk:10.0f} {s["dc"]:10.4f} {s["onoff"]:10.4f} {s["fm"]:10.4f}')
