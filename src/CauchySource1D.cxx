#include "CauchySource1D.hxx"

namespace JJCorrFitter
{
    CauchySource1D::CauchySource1D() : m_invariantRadius(2.f)
    {
        m_sourceFunctionName = "Cauchy";
        m_numberOfParams = 1;
    }

    void CauchySource1D::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("CauchySource1D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_invariantRadius = pars.at(0);
    }

    double CauchySource1D::GetValue(float rStar) const noexcept
    {
        return exp(-1. * rStar / (m_invariantRadius/*  * m_invariantRadius */));
    }
} // namespace JJCorrFitter
