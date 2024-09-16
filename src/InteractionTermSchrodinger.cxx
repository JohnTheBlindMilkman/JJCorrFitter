#include "InteractionTermSchrodinger.hxx"

namespace JJCorrFitter
{
    InteractionTermSchrodinger::InteractionTermSchrodinger() : 
    m_waveFunction(new CWaveFunction_pp_schrod("./wfparameters.dat"))
    {
        m_nqMax = m_waveFunction->GetNQMAX();
    }

    InteractionTermSchrodinger::~InteractionTermSchrodinger()
    {}

    InteractionTermSchrodinger::InteractionTermSchrodinger(InteractionTermSchrodinger &&other) noexcept
    {
        m_waveFunction = std::move(other.m_waveFunction);
        other.m_waveFunction = nullptr;
        m_nqMax = std::move(other.m_nqMax);
    }

    InteractionTermSchrodinger& InteractionTermSchrodinger::operator=(InteractionTermSchrodinger &&other) noexcept
    {
        m_waveFunction = std::move(other.m_waveFunction);
        other.m_waveFunction = nullptr;
        m_nqMax = std::move(other.m_nqMax);

        return *this;
    }

    void InteractionTermSchrodinger::SetParameters(float kStar, float cosTheta) noexcept
    {
        m_kStar = kStar;
        m_cosTheta = cosTheta;
    }

    double InteractionTermSchrodinger::GetValue(float rStar)
    {
        return m_waveFunction->GetPsiSquared(m_kStar,rStar,m_cosTheta);
    }

} // namespace JJCorrFitter
