import ROOT
from array import array

w0 = ROOT.TMath.TwoPi()

def WaveletFD(x, x0, s):
    return ROOT.TMath.Exp(-1 * ROOT.TMath.Power((x - x0) / s, 2) / 2)


def CWT_mw(FileName):
    DFtime = array('d')
    DFflux = array('d')

    DataFile = open(FileName, 'r')
    for line in DataFile:
        if line[0] != '#':
            try:
                values = [float(val) for val in line.split()]
                tmp = values[0] / 86400 + 51910.5
                DFtime.append(tmp)
                DFflux.append(values[2])
            except Exception as e:
                print(f"Error reading data file: {e}")
                return

    n = len(DFflux)
    zeros = array('d', [0.0] * n)
    Period = DFtime[-1] - DFtime[0]

    FMed = ROOT.TGraph(len(DFtime), array('d', DFtime), array('d', DFflux))

    DFfluxFFTDataRE = array('d', [0.0] * n)
    DFfluxFFTDataIM = array('d', [0.0] * n)
    DFfreqFFTData = array('d', [0.0] * n)
    mfftre = array('d', [0.0] * n)
    CWTfftre = array('d', [0.0] * n)
    CWTfftim = array('d', [0.0] * n)
    CWTre = array('d', [0.0] * n)
    CWTim = array('d', [0.0] * n)

    fft_data = ROOT.TFFTRealComplex(n)
    if not fft_data:
        return

    fft_data.SetPoints(DFflux.buffer_info()[0])
    fft_data.Transform()
    fft_data.GetPointsComplex(DFfluxFFTDataRE, DFfluxFFTDataIM)
    for i in range(n):
        DFfreqFFTData[i] = i

    TMin = 2 * Period / n
    TMax = Period
    TStart = 90
    TStop = 1500
    TNumber = 200
    T = TMin
    Tbins = array('d', [0.0] * (TNumber + 1))
    scaleFactor = 2 * TNumber / n
    TCbins = array('d', [0.0] * TNumber)

    FStart = 1 / TStop
    FStop = 1 / TStart
    FStep = (FStop - FStart) / TNumber
    F = array('d', [0.0] * TNumber)

    for j in range(TNumber):
        F[j] = FStart + j * FStep
        TCbins[TNumber - j - 1] = 1 / F[j]
        Tbins[TNumber - j] = 1 / (FStart + (2 * j - 1) * FStep / 2)
    Tbins[0] = 1 / (FStart + (2 * TNumber - 1) * FStep / 2)

    Proy = array('d', [0.0] * TNumber)

    aHist2 = ROOT.TH2D("hcwt", "CWT", n, DFtime[0], DFtime[-1], TNumber, Tbins)

    for j in range(1, TNumber + 1):
        for i in range(n):
            mfftre[i] = WaveletFD(i, Period * F[j], 1)

        for i in range(n):
            CWTfftre[i] = DFfluxFFTDataRE[i] * mfftre[i]
            CWTfftim[i] = DFfluxFFTDataIM[i] * mfftre[i]

        fft_CWT = ROOT.TVirtualFFT.FFT(1, array('I', [n]), "C2CBACKWARD P")
        if not fft_CWT:
            return

        fft_CWT.SetPointsComplex(CWTfftre, CWTfftim)
        fft_CWT.Transform()

        fft_CWT.GetPointsComplex(CWTre, CWTim)

        Proy[TNumber - j] = 0

        for i in range(n):
            val = ROOT.TMath.Sqrt(TCbins[TNumber - j]) * ROOT.TMath.Sqrt(CWTre[i] * CWTre[i] + CWTim[i] * CWTim[i])
            aHist2.SetBinContent(i + 1, TNumber - j + 1, val)
            Proy[TNumber - j] += val
        Proy[TNumber - j] /= n

    MyCanvas = ROOT.TCanvas("C1", FileName, 720, 1080)
    MyCanvas.Divide(1, 3)
    MyCanvas.cd(1)
    FMed.SetTitle(FileName)
    FMed.SetMarkerSize(1.5)
    FMed.SetMarkerStyle(21)
    FMed.Draw("APL")
    FMed.Draw()
    MyCanvas.cd(2)
    aHist2.SetStats(ROOT.kFALSE)
    ROOT.gPad.SetLogy()
    aHist2.GetYaxis().SetMoreLogLabels()
    aHist2.GetYaxis().SetNoExponent()
    aHist2.Draw("colz")
    MyCanvas.cd(3)

    PGraph = ROOT.TGraph(TNumber, TCbins, Proy)
    PGraph.SetTitle("CWT average over measured period")
    ROOT.gPad.SetLogx()
    PGraph.GetXaxis().SetMoreLogLabels()
    PGraph.GetXaxis().SetNoExponent()
    PGraph.Draw()

    # Calculate the maximum period
    maxProy = 0.0
    keep_index = 0
    for i in range(TNumber):
        if Proy[i] > maxProy and TCbins[i] > 100:
            maxProy = Proy[i]
            keep_index = i
    print("Maximum period:", TCbins[keep_index])


    del DFfluxFFTDataRE
    del DFfluxFFTDataIM
    del DFfreqFFTData
    del mfftre
    del CWTfftre
    del CWTfftim
    del CWTre
    del CWTim
    del Tbins
    del TCbins

    del Proy


# Example usage:
CWT_mw("PKS0447-439/PKS0447-439_results_45dias_sorted_normalized.dat")
