#include "CorrelationFunction1D.hxx"

namespace JJCorrFitter
{
    /* CorrelationFunction1D::CorrelationFunction1D(const std::string &name, const std::string &title, float kStarMin, float kStarMax, int nPoints) : 
    m_ppSchroed(), m_source1D(), m_minKStar(kStarMin), m_maxKStar(kStarMax), m_nPoints(nPoints), 
    m_histogramName(name), m_histogramTitle(title), m_correlationFunctionPoints(0,nPoints)
    {
        m_numberOfParams = 2 + m_sourceFunction->GetNParams() + m_interactionTerm->GetNParams();
    } */

    CorrelationFunction1D::CorrelationFunction1D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact) :
    m_minKStar(0), m_maxKStar(500), m_nPoints(100), 
    m_histogramName("hCF_fit"), m_histogramTitle(""), m_correlationFunctionPoints(0,100)
    {
        m_sourceFunction = std::exchange(source,nullptr);
        m_interactionTerm = std::exchange(interact,nullptr);

        m_numberOfParams = 2;
        m_totalNumberOfParams = m_numberOfParams + m_sourceFunction->GetNParams() + m_interactionTerm->GetNParams();
    }

    double CorrelationFunction1D::CalculatePoint()
    {
        std::cerr << "CorrelationFunction1D::CalculatePoint - This is a deleted function" << std::endl; 
        return 0;
    }

    double CorrelationFunction1D::CalculatePoint(float kStar)
    {
        m_interactionTerm->SetMomentum(kStar);

        // implemented as in kernel.cc for identical particles in CorAL, i.e. Koonin-Pratt, but instead of |Psi|^2 we use Kernel in 1D
        ROOT::Math::WrappedMultiFunction wfKP([&](const double *x)
        {
            return m_sourceFunction->GetValue(x[0]) * m_interactionTerm->GetValue(x[0],x[1]) * (gsl_sf_legendre_Pl(0,x[1]) + gsl_sf_legendre_Pl(2,x[1])); 
        },2);
        ROOT::Math::IntegratorMultiDim integKP(wfKP,ROOT::Math::IntegrationMultiDim::kDEFAULT,-1.0,-1.0,100000);


        const std::vector<double> integMin = {std::numeric_limits<float>::epsilon(),-1};
        const std::vector<double> integMax = {500,1};
        return integKP.Integral(integMin.data(),integMax.data());
    }

    std::unique_ptr<TH1> CorrelationFunction1D::Evaluate()
    {
        std::unique_ptr<TH1D> hist = std::make_unique<TH1D>(m_histogramName.data(),m_histogramTitle.data(),m_nPoints,m_minKStar,m_maxKStar);

        for (int bin = 1; bin <= m_nPoints; ++bin)
            hist->SetBinContent(bin,m_corrFuncParams.at(0) * (m_corrFuncParams.at(1) * (CalculatePoint(hist->GetBinCenter(bin)) - 1) + 1));

        return hist;
    }

    void CorrelationFunction1D::SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars)
    {
        if (generalPars.size() != m_numberOfParams)
        {
            throw std::length_error("CorrelationFunction1D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_corrFuncParams = generalPars;
        m_sourceFunction->SetParameters(srcPars);
        m_interactionTerm->SetParameters(psiPars);
    }
} // namespace JJCorrFitter
