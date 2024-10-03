#include "SourceFunction1D.hxx"

namespace JJCorrFitter
{
    SourceFunction1D::SourceFunction1D() : m_invariantRadius(2.f)
    {
        m_sourceFunctionName = "Gaussian";
        m_numberOfParams = 1;
    }

    void SourceFunction1D::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("SourceFunction1D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_invariantRadius = pars.at(0);
    }

    double SourceFunction1D::GetValue(float rStar) const noexcept
    {
        return pow(4 * ROOT::Math::Pi() * m_invariantRadius * m_invariantRadius,-1.5) * exp(-0.25*rStar*rStar/(m_invariantRadius*m_invariantRadius));
    }
} // namespace JJCorrFitter
