#include "InteractionTermSchrodinger.hxx"

namespace JJCorrFitter
{
    InteractionTermSchrodinger::InteractionTermSchrodinger() : 
    m_waveFunction(new CWaveFunction_pp_schrod(m_coralParFile.data()))
    {
        m_nqMax = m_waveFunction->GetNQMAX();
    }

    InteractionTermSchrodinger::~InteractionTermSchrodinger()
    {}

    InteractionTermSchrodinger::InteractionTermSchrodinger(InteractionTermSchrodinger &&other) noexcept
    {
        other.m_waveFunction = std::move(m_waveFunction);
        m_waveFunction = nullptr;
        other.m_nqMax = std::move(m_nqMax);
    }

    InteractionTermSchrodinger& InteractionTermSchrodinger::operator=(InteractionTermSchrodinger &&other) noexcept
    {
        other.m_waveFunction = std::move(m_waveFunction);
        m_waveFunction = nullptr;
        other.m_nqMax = std::move(m_nqMax);
    }

    double InteractionTermSchrodinger::GetValue(int kStar, float rStar, float cosTheta)
    {
        if (kStar < 0 || kStar > m_nqMax)
        {
            throw std::runtime_error("Provided value of kStar is outside of the bounds given by the CWaveFunction_pp_schrod object");
        }
        return m_waveFunction->CalcPsiSquared(kStar,rStar,cosTheta);
    }

} // namespace JJCorrFitter
