#include "InteractionTermSchrodinger.hxx"

namespace JJCorrFitter
{
    InteractionTermSchrodinger::InteractionTermSchrodinger() : 
        m_waveFunction("/home/jedkol/Downloads/JJCorrFitter/build/tests/wfparameters.dat"), m_kStar(5.)
    {
        m_InteractionTermName = "p-p Schroedinger solution (CorAL)";
        m_numberOfParams = 0;
        m_nqMax = m_waveFunction.GetNQMAX();
    }

    InteractionTermSchrodinger::~InteractionTermSchrodinger()
    {}

    InteractionTermSchrodinger::InteractionTermSchrodinger(InteractionTermSchrodinger &&other) noexcept : 
        m_waveFunction(std::move(other.m_waveFunction)), m_nqMax(std::move(other.m_nqMax))
    {}

    InteractionTermSchrodinger& InteractionTermSchrodinger::operator=(InteractionTermSchrodinger &&other) noexcept
    {
        m_waveFunction = std::move(other.m_waveFunction);
        m_nqMax = std::move(other.m_nqMax);

        return *this;
    }

    void InteractionTermSchrodinger::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("InteractionTermSchrodinger::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }
    }

    void InteractionTermSchrodinger::SetMomentum(float kStar)
    {
        m_kStar = kStar;
    }

    void InteractionTermSchrodinger::SetNMomentumBins(std::vector<float> nBins)
    {
        m_nMomBins = nBins;
        m_grid = Grid<double>(m_nMomBins.size(),m_nRPoints,m_nCosThetaPoints);
        PopulateGrid();
    }
    
    double InteractionTermSchrodinger::GetValue(float rStar, float cosTheta)
    {
        std::size_t momBin = std::distance(m_nMomBins.begin(),std::lower_bound(m_nMomBins.begin(),m_nMomBins.end(),m_kStar));
        std::size_t rBin = std::distance(m_rPoints.begin(),std::lower_bound(m_rPoints.begin(),m_rPoints.end(),rStar));
        std::size_t ctBin = std::distance(m_cosThetaPoints.begin(),std::lower_bound(m_cosThetaPoints.begin(),m_cosThetaPoints.end(),cosTheta));

        return m_grid.at(momBin,rBin,ctBin);
    }

    void InteractionTermSchrodinger::PopulateGrid()
    {
        for (std::size_t qBin = 0; qBin < m_nMomBins.size(); ++qBin)
            for (std::size_t rBin = 0; rBin < m_nRPoints; ++rBin)
                for (std::size_t ctBin = 0; ctBin < m_nCosThetaPoints; ++ctBin)
                    m_grid(qBin,rBin,ctBin) = m_waveFunction.GetPsiSquared(m_nMomBins.at(qBin),m_rPoints.at(rBin),m_cosThetaPoints.at(ctBin));
    }

} // namespace JJCorrFitter
