#include "InteractionTermSchrodinger.hxx"

namespace JJCorrFitter
{
    InteractionTermSchrodinger::InteractionTermSchrodinger() : 
    m_waveFunction(new CWaveFunction_pp_schrod(m_coralParFile.data())), 
    m_nqMax(m_waveFunction->GetNQMAX())
    {
    }

    InteractionTermSchrodinger::~InteractionTermSchrodinger()
    {}

    double InteractionTermSchrodinger::GetValue(int kStar, float rStar, float cosTheta)
    {
        if (kStar < 0 || kStar > m_nqMax)
        {
            std::cerr << "InteractionTermSchrodinger::GetValue - kStar is out of bounds (0,NQMAX)" << std::endl;
            std::exit(1);
        }
        return m_waveFunction->CalcPsiSquared(kStar,rStar,cosTheta);
    }

} // namespace JJCorrFitter
