#include "CorrelationFunction1D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction1D::CorrelationFunction1D(float kStarMin, float kStarMax, int nPoints) : 
    m_ppSchroed(), m_source1D(), m_cosTheta(1), m_minKStar(kStarMin), m_maxKStar(kStarMax), 
    m_nPoints(nPoints), m_correlationFunctionPoints(0,nPoints)
    {}

    double CorrelationFunction1D::CalculatePoint()
    {
        std::cerr << "CorrelationFunction1D::CalculatePoint - This is a deleted function" << std::endl; 
        return 0;
    }

    double CorrelationFunction1D::CalculatePoint(float kStar)
    {
        m_ppSchroed.SetParameters(kStar,m_cosTheta);

        ROOT::Math::WrappedFunction wf([&](double x){return m_source1D.GetValue(x) * m_ppSchroed.GetValue(x); });
        ROOT::Math::Integrator integ(wf,ROOT::Math::IntegrationOneDim::kGAUSS,-1.0,-1.0,1000);

        return integ.IntegralUp(std::numeric_limits<float>::epsilon());
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
