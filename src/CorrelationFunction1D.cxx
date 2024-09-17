#include "CorrelationFunction1D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction1D::CorrelationFunction1D(const std::string &name, const std::string &title, float kStarMin, float kStarMax, int nPoints) : 
    m_ppSchroed(), m_source1D(), m_minKStar(kStarMin), m_maxKStar(kStarMax), m_nPoints(nPoints), 
    m_histogramName(name), m_histogramTitle(title), m_correlationFunctionPoints(0,nPoints)
    {}

    double CorrelationFunction1D::CalculatePoint()
    {
        std::cerr << "CorrelationFunction1D::CalculatePoint - This is a deleted function" << std::endl; 
        return 0;
    }

    double CorrelationFunction1D::CalculatePoint(float kStar)
    {
        m_ppSchroed.SetParameters(kStar);

        // implemented as in kernel.cc for identical particles in CorAL, i.e. Koonin-Pratt, but instead of |Psi|^2 we use Kernel in 1D
        ROOT::Math::WrappedMultiFunction wfKP([&](double *x)
        {
            return m_source1D.GetValue(x[0]) * (m_ppSchroed.GetValue(x[0],x[1]) - 1) * (gsl_sf_legendre_Pl(0,x[1]) + gsl_sf_legendre_Pl(2,x[1])); 
        });
        ROOT::Math::IntegratorMultiDim integKP(wfKP,ROOT::Math::IntegrationMultiDim::kDEFAULT,-1.0,-1.0,1000);


        const std::vector<double> integMin = {std::numeric_limits<float>::epsilon(),-1};
        const std::vector<double> integMax = {500,1};
        return integKP.Integral(integMin.data(),integMax.data());
    }

    std::unique_ptr<TH1> CorrelationFunction1D::Evaluate()
    {
        std::unique_ptr<TH1D> hist = std::make_unique<TH1D>(m_histogramName.data(),m_histogramTitle.data(),m_nPoints,m_minKStar,m_maxKStar);

        for (int bin = 1; bin <= m_nPoints; ++bin)
            hist->SetBinContent(bin,CalculatePoint(hist->GetBinCenter(bin)));

        return hist;
    }

    void CorrelationFunction1D::SetIntegrationRange(float rStarMin, float rStarMax) noexcept
    {
        m_minRStar = rStarMin;
        m_MaxRStar = rStarMax;
    }

    void CorrelationFunction1D::SetParameters() noexcept
    {
        std::cerr << "CorrelationFunction1D::CalculatePoint - This is a deleted function" << std::endl; 
    }

    void CorrelationFunction1D::SetParameters(float rInv) noexcept
    {
        m_source1D.SetParameters(rInv);
    }
} // namespace JJCorrFitter
