#include "TFile.h"
#include "CorrelationFunction1D.hxx"

int main()
{
    JJCorrFitter::CorrelationFunction1D cf("hCF1D","One-dimensional correlation function",5.,200.,100);
    cf.SetParameters(5.);

    TFile *otp = TFile::Open("test.root","recreate");

    cf.Evaluate()->Write();

    otp->Close();
}