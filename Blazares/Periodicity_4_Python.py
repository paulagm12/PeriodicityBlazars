import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.signal import morlet
from scipy.constants import pi
import pywt

w0 = 2 * np.pi

def WaveletFD(x, x0, s):
    return np.exp(-1 * ((x - x0) / s) ** 2 / 2)

def CWT_mw(file_name):
    # Carga de datos desde el archivo
    DFtime, sDFtime, DFflux, sDFflux = [], [], [], []

    with open(file_name, 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                try:
                    c0, c1, c2, c3, c4, c5, c6, c7 = map(float, line.split())
                    tmp = c0 / 86400 + 51910.5  # Converts MET to Julian days
                    DFtime.append(tmp)
                    DFflux.append(c2)
                except Exception as e:
                    print("Error reading data file:", str(e))
                    return
            else:
                print("Skipping line")

    # Display loaded data size and range
    print(f'DFtime size: {len(DFtime)}, Start: {DFtime[0]}, End: {DFtime[-1]}')

    n = len(DFflux)
    zeros = np.zeros(n)
    Period = DFtime[-1] - DFtime[0]

    FMed = plt.figure()
    # Graphs light curve
    plt.plot(DFtime, DFflux, marker='s', linestyle='-', color='black')
    plt.xlabel('Time (JD)')
    plt.ylabel('Flux (ph cm-2 s-1)')
    plt.title(file_name)
    plt.show() # This will have to be eliminated and put into a subplot
# HASTA AQUÍ TODO BIEN TODO LO DE ARRIBA BIEN

    # FFT
    DFfluxFFT = fft(DFflux)
    DFfluxFFTDataRE = DFfluxFFT.real
    DFfluxFFTDataIM = DFfluxFFT.imag

    # Setting up wavelet transform parameters
    TMin = 2 * Period / n
    TMax = Period
    TStart = 90
    TStop = 1500
    TNumber = 200
    Tbins = np.zeros(TNumber + 1)
    FStart = 1 / TStop
    FStop = 1 / TStart
    FStep = (FStop - FStart) / TNumber
    F = np.zeros(TNumber)

    # Frequency and period bins
    for j in range(TNumber):
        F[j] = FStart + j * FStep
        Tbins[TNumber - j] = 1 / (FStart + (2 * j - 1) * FStep / 2)
    Tbins[0] = 1 / (FStart + (2 * TNumber - 1) * FStep / 2)

    # Placeholder for results
    CWT_re = np.zeros((TNumber, n))
    CWT_im = np.zeros((TNumber, n))

    # Perform wavelet transform (simplified)
    for j in range(1, TNumber + 1):
        mfftre = WaveletFD(np.arange(n), Period * F[j - 1], 1)
        CWTfftre = DFfluxFFTDataRE * mfftre
        CWTfftim = DFfluxFFTDataIM * mfftre

        # Inverse FFT
        CWT = ifft(CWTfftre + 1j * CWTfftim)
        CWT_re[j - 1, :] = CWT.real
        CWT_im[j - 1, :] = CWT.imag

    # Visualization
    fig, ax = plt.subplots()
    im = ax.imshow(np.sqrt(CWT_re**2 + CWT_im**2), aspect='auto', extent=[DFtime[0], DFtime[-1], Tbins[-1], Tbins[0]], cmap='jet')
    #ax.set_yscale('log')
    plt.colorbar(im, ax=ax)
    plt.show()


CWT_mw("/home/paula/Documents/Internship/AstroECFM/Blazares/PKS0447-439/PKS0447-439_results_45dias_sorted_normalized.dat")
