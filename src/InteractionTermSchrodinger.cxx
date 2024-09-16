#include "InteractionTermSchrodinger.hxx"

namespace JJCorrFitter
{
    InteractionTermSchrodinger::InteractionTermSchrodinger() : 
    m_waveFunction(new CWaveFunction_pp_schrod(m_coralParFile.data()))
    {
        m_nqMax = m_waveFunction->GetNQMAX();
        m_kStar = 1;
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

    void InteractionTermSchrodinger::SetParameters(int kStar) noexcept
    {
        m_kStar = std::min(std::max(kStar,m_kStarMin),m_nqMax); // take the smaller value: tmp or nQmax (limited by CorAL); where tmp is greater value between k* and m_kStarMin
        // baisically takes kStar between m_kStarMin and nQmax or rounds up/down the set value of kStar to match the limits
    }

    double InteractionTermSchrodinger::GetValue(float rStar, float cosTheta)
    {
        return m_waveFunction->CalcPsiSquared(m_kStar,rStar,cosTheta);
    }

} // namespace JJCorrFitter
