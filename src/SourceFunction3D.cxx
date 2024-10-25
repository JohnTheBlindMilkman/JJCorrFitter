#include "SourceFunction3D.hxx"

namespace JJCorrFitter
{
    SourceFunction3D::SourceFunction3D() : m_outRadius(2.f), m_sideRadius(2.f), m_longRadius(2.f)
    {
        m_sourceFunctionName = "Gaussian";
        m_numberOfParams = 3;
    }

    void SourceFunction3D::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("SourceFunction3D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_outRadius = pars.at(0);
        m_sideRadius = pars.at(1);
        m_longRadius = pars.at(2);
    }

    double SourceFunction3D::GetValue(float rOut, float rSide, float rLong) const noexcept
    {
        return exp(-0.25 * ((rOut * rOut / (m_outRadius * m_outRadius)) + (rSide * rSide / (m_sideRadius * m_sideRadius)) + (rLong * rLong / (m_longRadius * m_longRadius))));
    }
} // namespace JJCorrFitter
