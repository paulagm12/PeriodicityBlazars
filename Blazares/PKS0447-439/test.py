import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

w0 = 2 * np.pi

def wavelet_fd(x, x0, s):
    return np.exp(-1 * ((x - x0) / s) ** 2 / 2)

def cwt_mw(file_name):
    with open(file_name, 'r') as data_file:
        DFtime, DFflux = [], []
        for line in data_file:
            if not line.startswith('#'):
                c0, c1, c2, c3, c4, c5, c6, c7 = map(float, line.split())
                tmp = c0 / 86400 + 51910.5
                DFtime.append(tmp)
                DFflux.append(c2)

    n = len(DFflux)
    zeros = np.zeros(n)
    Period = DFtime[-1] - DFtime[0]

    FMed = plt.figure()
    plt.plot(DFtime, DFflux, marker='o', linestyle='-')
    plt.title(file_name)
    plt.show()

    DFflux_fft = fft(DFflux)
    DFfreqFFTData = np.fft.fftfreq(n)

    TMin = 2 * Period / n
    TMax = Period
    TStart = 90
    TStop = 1500
    TNumber = 200
    T = TMin
    Tbins = np.zeros(TNumber + 1)
    TCbins = np.zeros(TNumber)

    FStart = 1 / TStop
    FStop = 1 / TStart
    FStep = (FStop - FStart) / TNumber
    F = np.zeros(TNumber)

    for j in range(TNumber):
        F[j] = FStart + j * FStep
        TCbins[TNumber - j - 1] = 1 / F[j]
        Tbins[TNumber - j] = 1 / (FStart + (2 * j - 1) * FStep / 2)
    Tbins[0] = 1 / (FStart + (2 * TNumber - 1) * FStep / 2)

    Proy = np.zeros(TNumber)

    for j in range(TNumber):
        mfftre = wavelet_fd(np.arange(n), Period * F[j], 1)

        CWTfftre = DFflux_fft * mfftre

        CWT = ifft(CWTfftre)

        CWTre = np.real(CWT)
        CWTim = np.imag(CWT)

        for i in range(n):
            Proy[TNumber - j - 1] += np.sqrt(TCbins[TNumber - j - 1]) * np.sqrt(CWTre[i] ** 2 + CWTim[i] ** 2)

    Proy /= n

    CWT = np.abs(CWT)  # Take the absolute values of CWT if not already done

    # Reshape CWT to 2D array for visualization
    CWT = np.tile(CWT, (len(TCbins), 1)).T  # Make the data 2D for imshow

    # Plot the CWT
    plt.figure()
    plt.imshow(CWT, aspect='auto', origin='lower', extent=[DFtime[0], DFtime[-1], FStart, FStop])
    plt.colorbar()
    plt.title('CWT')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.show()

    # Plot CWT average over measured period
    plt.figure()
    plt.plot(TCbins, Proy)
    plt.xscale('log')
    plt.title('CWT average over measured period')
    plt.show()

    # Calculate the maximum period
    max_proy = 0.0
    keep_index = 0
    for i in range(TNumber):
        if Proy[i] > max_proy and TCbins[i] > 100:
            max_proy = Proy[i]
            keep_index = i
    print(f"Maximum period: {TCbins[keep_index]}")

cwt_mw("PKS0447-439_results_45dias_sorted_normalized.dat")
# cwt_mw("OJ287/OJ287_normalized.dat")
