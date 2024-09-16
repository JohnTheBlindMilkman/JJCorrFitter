#include "SourceFunction1D.hxx"

namespace JJCorrFitter
{
    SourceFunction1D::SourceFunction1D() : m_invariantRadius(2.f)
    {}

    void SourceFunction1D::SetParameters(float rInv) noexcept
    {
        m_invariantRadius = rInv;
    }

    double SourceFunction1D::GetValue(float rStar) const noexcept
    {
        return pow(4 * ROOT::Math::Pi() * m_invariantRadius * m_invariantRadius,-1.5) * exp(-0.25*rStar*rStar/(m_invariantRadius*m_invariantRadius));
    }
} // namespace JJCorrFitter
