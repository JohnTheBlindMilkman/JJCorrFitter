#include "SourceFunction1D.hxx"

namespace JJCorrFitter
{
    double SourceFunction1D::GetValue(float rStar, float rInv) const noexcept
    {
        return pow(4 * ROOT::Math::Pi() * rInv * rInv,-1.5) * exp(-0.25*rStar*rStar/(rInv*rInv));
    }
} // namespace JJCorrFitter
