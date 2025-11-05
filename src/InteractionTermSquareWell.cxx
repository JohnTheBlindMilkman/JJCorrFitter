#include "InteractionTermSquareWell.hxx"

namespace JJCorrFitter
{
    InteractionTermSquareWell::InteractionTermSquareWell() : m_kStar(5.)
    {
        m_InteractionTermName = "p-p square well";
        m_numberOfParams = 0;
    }

    void InteractionTermSquareWell::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("InteractionTermSquareWell::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }
    }

    void InteractionTermSquareWell::SetMomentum(double kStar) noexcept
    {
        m_kStar = kStar;
    }
    
    double InteractionTermSquareWell::GetValue(double rStar,double cosTheta) noexcept
    {
        return 1 - 0.5 * std::cos(2 * rStar * m_kStar * cosTheta / m_hBarC);
    }

} // namespace JJCorrFitter
