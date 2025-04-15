#include "DoubleGaussian1D.hxx"

namespace JJCorrFitter
{
    DoubleGaussian1D::DoubleGaussian1D() : m_invariantRadius1(2.f), m_invariantRadius2(2.f)
    {
        m_sourceFunctionName = "Double Gaussian";
        m_numberOfParams = 2;
    }

    void DoubleGaussian1D::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("DoubleGaussian1D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_invariantRadius1 = pars.at(0);
        m_invariantRadius2 = pars.at(1);
    }

    double DoubleGaussian1D::GetValue(float rStar) const noexcept
    {
        return exp(-0.5 * rStar * rStar / (m_invariantRadius1 * m_invariantRadius1)) + exp(-0.5 * rStar * rStar / (m_invariantRadius2 * m_invariantRadius2));
    }
} // namespace JJCorrFitter
