import numpy as np
import cmath


def DFT(t, x_t):
    """
    Computes the Discrete Fourier Transform (DFT) of a signal.

    Parameters:
        t (array): time samples (uniformly spaced)
        x_t (array): signal samples
    Returns:
        f (array): frequency bins (Hz)
        X_f (array): DFT coefficients
    """
    N = np.shape(t)[0]  # number of samples
    fs = 1 / (t[1] - t[0])  # sampling frequency

    X_f = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            X_f[k] += x_t[n] * cmath.exp(complex(0, -2*np.pi/N * k*n))
    f = np.zeros(N)

    # Only take positive frequencies
    half_N = N // 2
    f = np.arange(half_N) * fs / N
    X_f = X_f[:half_N]

    return f, X_f


def IDFT(f, X_f):
    # optional
    """
    Computes the Inverse Discrete Fourier Transform (IDFT) of a signal.

    Parameters:
        f (array): frequency bins (Hz)
        X_f (array): DFT coefficients
    Returns:
        t (array): time samples (uniformly spaced)
        x_t (array): signal samples
    """
    N = np.shape(f)[0]
    Ts = np.abs(1 / (N*(f[1] - f[0])))
    x_t = np.zeros(N, dtype=complex)
    t = np.zeros(N)

    for n in range(N):
        t[n] = n*Ts
        for k in range(N):
            x_t[n] += X_f[k] * cmath.exp(complex(0, 2 * np.pi / N * k * n))
        x_t[n] = x_t[n] / N

    return t, np.real(x_t)
