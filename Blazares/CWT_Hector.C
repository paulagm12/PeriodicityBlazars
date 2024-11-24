// root script file

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TMultiGraph.h"
#include "TH2D.h"
#include "TVirtualFFT.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>


const Double_t w0 = TMath::TwoPi();
//Double_t max_prob = 0.0;

Double_t WaveletFD(Double_t x, Double_t x0, Double_t s)
{
    return TMath::Exp( -1*TMath::Power( (x-x0)/s , 2 )/2 ) ;
}


void CWT_mw(const char * FileName)
{
    /* Carga de datos desde el archivo */
    std::vector<Double_t> DFtime, sDFtime, DFflux, sDFflux;
    std::ifstream DataFile(FileName);
    Int_t lineNum = 0;
    Double_t c0, c1, c2, c3, c4, c5, c6, c7,tmp;
    char aBuffer[512];
    while( !DataFile.eof() )
    {
        
        DataFile.getline(aBuffer,512);

        // std::cout << aBuffer << std::endl;

        if( aBuffer[0] != '#' )
        {
            std::stringstream input_stringstream(aBuffer);
            try {
                input_stringstream >> c0 >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7;
                tmp = c0/86400 + 51910.5; //Convirtiendo de MET a algo juliano...
                DFtime.push_back( tmp );
                DFflux.push_back( c2 );
            }
            catch (...){
                std::cout << "Error reading data file" << std::endl;
                return;
            }
        }
        else
            std::cout << "skipping line " << std::endl;
        
        
    }

    std::cout << DFtime.size() << " " << DFtime.back() << " " << DFtime.front() << std::endl;


    

    Int_t n = DFflux.size(); //Número de récords
    Double_t *zeros = new Double_t[n];
    Double_t Period = DFtime.back() - DFtime.front(); // Calcular el tiempo total de estudio
    for (Int_t i=0; i<n; i++)
    {
        zeros[i] = 0;
    }

    // Int_t n = 1024;
    // Double_t Period = 100;

    // Double_t *zeros = new Double_t[n];

    // DFflux.clear();
    // DFtime.clear();
    // for (Int_t i=0; i<n; i++)
    // {
    //     zeros[i] = 0;
    //     tmp = double(i)/(n+1)*Period;
    //     DFflux.push_back( ( tmp > 50 ? TMath::Cos(TMath::TwoPi()*(tmp/2)) : TMath::Cos(TMath::TwoPi()*(tmp/10)) )  ) ;
    //     // DFflux.push_back( ( tmp > 50 ? 1 : 0  ) );
    //     // DFflux.push_back( TMath::Cos(TMath::TwoPi()*(tmp/50)) + 0*TMath::Cos(TMath::TwoPi()*(tmp/10))   ) ;
    //     DFtime.push_back(tmp);
    // }

    TGraph *FMed = new TGraph(DFtime.size(),DFtime.data(),DFflux.data()); // Graficando curva de luz

    Double_t *DFfluxFFTDataRE = new Double_t[n];
    Double_t *DFfluxFFTDataIM = new Double_t[n];
    Double_t *DFfreqFFTData = new Double_t[n];
    Double_t *mfftre = new Double_t[n];
    Double_t *CWTfftre = new Double_t[n];
    Double_t *CWTfftim = new Double_t[n];
    Double_t *CWTre = new Double_t[n];
    Double_t *CWTim = new Double_t[n];

    TVirtualFFT *fft_data = TVirtualFFT::FFT(1, &n, "C2CFORWARD P");
    if (!fft_data) return;
    fft_data->SetPointsComplex(DFflux.data(),zeros);
    fft_data->Transform();

    for (Int_t i=0; i<n; i++){
        fft_data->GetPointComplex(i, DFfluxFFTDataRE[i], DFfluxFFTDataIM[i]);
        DFfreqFFTData[i] =  i;
    }

    

    Double_t TMin = 2*Period/n;
    Double_t TMax = Period;
    Double_t TStart = 90;
    Double_t TStop = 1500;
    Int_t TNumber = 200;
    Double_t T = TMin;
    Double_t *Tbins = new Double_t[TNumber+1];
    Double_t scaleFactor = 2*TNumber/n;
    Double_t *TCbins = new Double_t[TNumber];

    Double_t FStart = 1/TStop;
    Double_t FStop = 1/TStart;
    Double_t FStep = ( FStop - FStart ) / TNumber;
    Double_t *F = new Double_t[TNumber];

    for( Int_t j=0; j<TNumber ; j++ )
    {
        F[j] = FStart + j*FStep;
        TCbins[TNumber-j-1] = 1 / F[j];
        Tbins[TNumber-j] = 1/ ( FStart + ( 2*j-1)*FStep/2 );
    }
    Tbins[0] = 1/ (FStart  + ( 2*TNumber-1)*FStep / 2 );

    Double_t *Proy = new Double_t[TNumber];

    TH2D *aHist2 = new TH2D("hcwt","CWT",n,DFtime.front(),DFtime.back(),TNumber,Tbins);

    
    for( Int_t j=1; j<=TNumber ; j++ )
    {

        for (Int_t i=0; i<n; i++){
            mfftre[i] = WaveletFD(i, Period*F[j], 1);
        }

        

        for (Int_t i=0; i<n; i++){
            CWTfftre[i] = DFfluxFFTDataRE[i]*mfftre[i] ;
            CWTfftim[i] = DFfluxFFTDataIM[i]*mfftre[i] ;
        }


        TVirtualFFT *fft_CWT = TVirtualFFT::FFT(1, &n, "C2CBACKWARD P");
        if (!fft_CWT) return;

        fft_CWT->SetPointsComplex(CWTfftre, CWTfftim);

        fft_CWT->Transform();

        fft_CWT->GetPointsComplex(CWTre,CWTim);
        
        Proy[TNumber - j ] = 0;
        
        for (Int_t i=0; i<n; i++)
        {
            aHist2->SetBinContent(i+1 , TNumber - j + 1,TMath::Sqrt(TCbins[TNumber - j]) * TMath::Sqrt(CWTre[i]*CWTre[i] + CWTim[i]*CWTim[i]));
            Proy[TNumber - j ] += TMath::Sqrt(TCbins[TNumber - j]) * TMath::Sqrt(CWTre[i]*CWTre[i] + CWTim[i]*CWTim[i]);
            std::cout << j << Proy[TNumber - j] << std::endl;
        }
        Proy[TNumber - j ] /= n;
        std::cout << "End" << std::endl;
    }
    
    

    TCanvas *MyCanvas = new TCanvas("C1",FileName,720,1080);
    MyCanvas->Divide(1,3);
    MyCanvas->cd(1);
    FMed->SetTitle(FileName);
    FMed->SetMarkerSize(1.5);
    FMed->SetMarkerStyle(21);
    FMed->Draw("APL");
    FMed->Draw();
    MyCanvas->cd(2);
    aHist2->SetStats(kFALSE);
    gPad->SetLogy();
    aHist2->GetYaxis()->SetMoreLogLabels();
    aHist2->GetYaxis()->SetNoExponent();
    aHist2->Draw("colz");
        

    MyCanvas->cd(3);

    TGraph *PGraph = new TGraph(TNumber,TCbins,Proy);
    PGraph->SetTitle("CWT average over measured period");
    gPad->SetLogx();
    PGraph->GetXaxis()->SetMoreLogLabels();
    PGraph->GetXaxis()->SetNoExponent();
    PGraph->Draw();
   


    delete [] DFfluxFFTDataRE;
    delete [] DFfluxFFTDataIM;
    delete [] DFfreqFFTData;
    delete [] mfftre;
    delete [] CWTfftre;
    delete [] CWTfftim;
    delete [] CWTre;
    delete [] CWTim;
    delete [] Tbins;
    delete [] TCbins;

    delete [] Proy;

}
