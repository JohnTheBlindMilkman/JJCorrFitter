#include "InteractionTermPhaseShift.hxx"

namespace JJCorrFitter
{
    InteractionTermPhaseShift::InteractionTermPhaseShift() : 
    m_waveFunction(new CWaveFunction_pp_phaseshift("./wfparameters.dat")), m_kStar(5.)
    {
        m_InteractionTermName = "p-p phase-shift (CorAL)";
        m_numberOfParams = 0;
        m_nqMax = m_waveFunction->GetNQMAX();
    }

    InteractionTermPhaseShift::~InteractionTermPhaseShift()
    {}

    InteractionTermPhaseShift::InteractionTermPhaseShift(InteractionTermPhaseShift &&other) noexcept
    {
        m_waveFunction = std::move(other.m_waveFunction);
        other.m_waveFunction = nullptr;
        m_nqMax = std::move(other.m_nqMax);
    }

    InteractionTermPhaseShift& InteractionTermPhaseShift::operator=(InteractionTermPhaseShift &&other) noexcept
    {
        m_waveFunction = std::move(other.m_waveFunction);
        other.m_waveFunction = nullptr;
        m_nqMax = std::move(other.m_nqMax);

        return *this;
    }

    void InteractionTermPhaseShift::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("InteractionTermPhaseShift::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }
    }

    void InteractionTermPhaseShift::SetMomentum(double kStar)
    {
        m_kStar = kStar;
    }
    
    double InteractionTermPhaseShift::GetValue(double rStar, double cosTheta)
    {
        return m_waveFunction->GetPsiSquared(m_kStar,rStar,cosTheta);
    }

} // namespace JJCorrFitter
