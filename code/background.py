"""
background.py
-------------
Excitation-wavelength dependence of the NON-SPIN-DEPENDENT background light
that lands inside the NV- detection window of a DAC experiment (Part E of
PLAN.md).

The sensitivity model in `nv_model.py` counts only photons emitted by NV
centers.  A real high-pressure ODMR setup collects, through the same
long-pass filter (typically >650 nm), light that has nothing to do with the
NV spin:

  * ruby R lines (694.2 / 692.8 nm) from the pressure calibrant -- these sit
    squarely INSIDE the NV detection window;
  * N3 / A-band luminescence of the anvil itself (nitrogen aggregates);
  * a broad, structureless anvil luminescence that grows with plastic
    deformation of the culet;
  * diamond Raman (1332 cm^-1) -- included only so that we can SHOW it does
    not reach the window for any of the excitation lines considered.

What matters for the green-vs-blue decision is not the absolute background
level (unknown without a measurement) but its EXCITATION-WAVELENGTH SHAPE.
Each channel k is therefore described by a normalised shape g_k(lambda_exc)
with g_k(532 nm) = 1, and the absolute level is carried by a single
dimensionless dial

    rho0 = B / R   evaluated at 532 nm excitation, ambient pressure

so the total background is

    B(lambda, P) = rho0 * R_ref * sum_k mix_k * g_k(lambda) * h_k(P)

with R_ref the model photon rate for green at ambient pressure.

Shape ingredients (all PHENOMENOLOGICAL, like the rate constants in
nv_model.py; the point of the analysis is that the CONCLUSION should not
depend on them in detail -- see the simplex randomisation in nv_bg.py):

  g_ruby : Cr3+ in Al2O3 has two broad absorption bands, U (~410 nm) and
           Y (~555 nm), with a well-known transmission window between them.
           470-500 nm therefore excites ruby POORLY -- a real advantage of
           blue excitation that the NV-only model cannot see.  Emission is
           the R line at ~694 nm, fully inside the detection window.
  g_n3   : N3 (ZPL 415 nm) absorbs increasingly strongly towards the UV and
           emits at 415-520 nm, i.e. almost entirely OUTSIDE a >650 nm
           window; only the weak red tail leaks in.  Penalises 405 nm.
  g_broad: deformation-related broad-band anvil luminescence.  Weak
           excitation-wavelength dependence, emission extending into the red
           so the window factor is large.  Optionally grows with pressure.
"""

import numpy as np

HBARC = 1239.84                    # eV * nm
LAM_REF = 532.0                    # normalisation wavelength for all shapes
WINDOW_NM = 650.0                  # long-pass edge of the NV detection window
RAMAN_SHIFT_CM = 1332.0            # first-order diamond Raman

CHANNELS = ('ruby', 'n3', 'broad')


# ---------------------------------------------------------------- shapes ---
def _gauss(lam, c, w):
    return np.exp(-0.5 * ((np.asarray(lam, float) - c) / w) ** 2)


def g_ruby_raw(lam):
    """Cr3+ absorption: U band ~410 nm + Y band ~555 nm, with a gap between.

    The gap (about 460-500 nm) is the physically interesting feature: blue
    excitation in that range barely pumps the ruby R line.
    """
    return 1.00 * _gauss(lam, 410.0, 33.0) + 1.15 * _gauss(lam, 555.0, 40.0) + 0.02


def g_n3_raw(lam):
    """N3 absorption rising towards the UV, times a small red-tail window factor.

    The absorption edge is modelled as an exponential (Urbach-like) rise
    below ~430 nm; the window factor is constant and small because N3 emits
    at 415-520 nm and only its tail reaches >650 nm.
    """
    lam = np.asarray(lam, float)
    absorb = np.exp(-(lam - 415.0) / 38.0)      # 1 at the ZPL, falls off to the red
    window_leak = 0.05                          # fraction of N3 emission above 650 nm
    return window_leak * absorb


def g_broad_raw(lam):
    """Deformation-related broad-band anvil luminescence.

    Mild monotonic increase towards the blue (shorter wavelength = more
    absorbing centres reachable), emission broad and extending into the red
    so essentially all of it is inside the window.
    """
    lam = np.asarray(lam, float)
    return 1.0 + 0.9 * (LAM_REF - lam) / LAM_REF


_RAW = {'ruby': g_ruby_raw, 'n3': g_n3_raw, 'broad': g_broad_raw}
_NORM = {k: float(f(LAM_REF)) for k, f in _RAW.items()}


def g_channel(name, lam):
    """Normalised background shape g_k(lambda_exc), with g_k(532 nm) = 1."""
    return _RAW[name](lam) / _NORM[name]


def default_mix():
    """Equal weights over the three channels (weights sum to 1)."""
    return {k: 1.0 / len(CHANNELS) for k in CHANNELS}


def pressure_factor(name, P, c_P=0.0):
    """Optional pressure growth of a channel (off by default).

    Only `broad` is given a pressure dependence: plastic deformation of the
    culet accumulates with load, so deformation luminescence is expected to
    grow.  `c_P` = relative increase at 120 GPa.
    """
    if name != 'broad' or c_P == 0.0:
        return np.ones_like(np.asarray(P, float))
    return 1.0 + c_P * np.clip(P, 0, None) / 120.0


def g_total(lam, P=0.0, mix=None, c_P=0.0):
    """Total normalised background shape summed over channels.

    Returns an array broadcast over `lam` and `P`.  Equals 1 at
    lambda = 532 nm, P = 0 for any mix whose weights sum to 1.
    """
    mix = default_mix() if mix is None else mix
    return sum(w * g_channel(k, lam) * pressure_factor(k, P, c_P)
               for k, w in mix.items() if w != 0.0)


# ----------------------------------------------------------- diamond Raman --
def raman_nm(lam_exc, order=1):
    """Stokes wavelength of order-n diamond Raman for a given excitation line."""
    lam_exc = np.asarray(lam_exc, float)
    return 1.0 / (1.0 / lam_exc - order * RAMAN_SHIFT_CM * 1e-7)


def raman_in_window(lam_exc, order=1, window=WINDOW_NM):
    """True if the order-n Raman line falls inside the detection window."""
    return raman_nm(lam_exc, order) >= window


def raman_report(lines=(405.0, 457.0, 473.0, 488.0, 532.0), window=WINDOW_NM):
    """Text table showing that diamond Raman never reaches the NV window."""
    rows = []
    for lam in lines:
        r1, r2 = float(raman_nm(lam, 1)), float(raman_nm(lam, 2))
        rows.append((lam, r1, r2, r1 >= window, r2 >= window))
    return rows


if __name__ == '__main__':
    print(f'detection window: > {WINDOW_NM:.0f} nm')
    print(f'{"exc":>6} {"1st Raman":>10} {"2nd Raman":>10}  in-window(1st/2nd)')
    for lam, r1, r2, i1, i2 in raman_report():
        print(f'{lam:6.0f} {r1:10.1f} {r2:10.1f}  {i1}/{i2}')
    print()
    for lam in (405, 457, 473, 475, 488, 532):
        parts = '  '.join(f'{k}={g_channel(k, lam):5.2f}' for k in CHANNELS)
        print(f'{lam:4d} nm  {parts}   total={float(g_total(lam)):5.2f}')
