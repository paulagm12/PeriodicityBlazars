import numpy as np
import matplotlib.pyplot as plt
import pywt
import pandas as pd

# Definir función de Transformada Wavelet Continua
def cwt(data, wavelet, scales):
    coefficients, frequencies = pywt.cwt(data, scales, wavelet)
    return coefficients, frequencies

# Cargar los datos desde el archivo
def load_data(file_path):
    # Leer el archivo ignorando líneas de comentarios
    data = pd.read_csv(file_path, sep='\s+', comment='#', header=None)
    return data

# Función principal para generar el scalogram
def plot_scalogram(file_path):
    # Cargar datos
    data = load_data(file_path)
    DFtime = data.iloc[:, 0].astype(float).values
    DFflux = data.iloc[:, 1].astype(float).values

    # Parámetros de la transformada wavelet continua
    TMin = 2 * (DFtime[-1] - DFtime[0]) / len(DFtime)
    TMax = DFtime[-1] - DFtime[0]
    TStart = 90
    TStop = 1500
    TNumber = 200
    scales = np.linspace(TMin, TMax, TNumber)

    # Realizar la transformada wavelet continua
    coefficients, frequencies = cwt(DFflux, 'morl', scales)

    # Crear el gráfico del scalogram
    plt.figure(figsize=(10, 6))

    # Gráfico de la señal original
    plt.subplot(2, 1, 1)
    plt.plot(DFtime, DFflux)
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.title('Original Signal')

    # Gráfico del scalogram
    plt.subplot(2, 1, 2)
    plt.imshow(np.abs(coefficients), extent=[DFtime[0], DFtime[-1], scales[-1], scales[0]], aspect='auto', cmap='jet')
    plt.colorbar(label='Magnitude')
    plt.xlabel('Time')
    plt.ylabel('Scale')
    plt.title('Scalogram')

    plt.tight_layout()
    plt.show()

# Ruta al archivo de datos
file_path = '/home/paula/Documents/Internship/AstroECFM/Blazares/PKS0447-439/PKS0447-439_results_45dias_sorted_normalized.dat'
plot_scalogram(file_path)
