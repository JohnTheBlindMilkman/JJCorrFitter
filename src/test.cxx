#include "TFile.h"
#include "CorrelationFunction1D.hxx"

int main()
{
    JJCorrFitter::CorrelationFunction1D cf(5.,200.,100);
    cf.SetParameters(1.);
    cf.SetIntegrationRange(0,500);

    TFile *otp = TFile::Open("test.root","recreate");

    cf.Evaluate()->Write();

    otp->Close();
    delete otp;
}