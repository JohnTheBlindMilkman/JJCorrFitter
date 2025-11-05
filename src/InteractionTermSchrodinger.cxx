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

    void InteractionTermSchrodinger::SetMomentum(double kStar)
    {
        m_kStar = kStar;
    }
    
    double InteractionTermSchrodinger::GetValue(double rStar, double cosTheta)
    {
        std::size_t momBin = std::distance(m_qPoints.begin(),std::lower_bound(m_qPoints.begin(),m_qPoints.end(),m_kStar));
        std::size_t rBin = std::distance(m_rPoints.begin(),std::lower_bound(m_rPoints.begin(),m_rPoints.end(),rStar));
        std::size_t ctBin = std::distance(m_ctPoints.begin(),std::lower_bound(m_ctPoints.begin(),m_ctPoints.end(),cosTheta));

        return m_grid.at(momBin,rBin,ctBin);
    }

    void InteractionTermSchrodinger::PopulateGrid()
    {
        auto qSize = m_qPoints.size();
        auto rSize = m_rPoints.size();
        auto ctSize = m_ctPoints.size();

        if (qSize == 0 || rSize == 0 || ctSize == 0)
            throw std::runtime_error("InteractionTermSchroedinger::PolulateGrid - current grid size is 0 in (at least) one dimension");

        m_grid = Grid<double>(qSize,rSize,ctSize);

        for (std::size_t qBin = 0; qBin < qSize; ++qBin)
            for (std::size_t rBin = 0; rBin < rSize; ++rBin)
                for (std::size_t ctBin = 0; ctBin < ctSize; ++ctBin)
                {
                    const double q = m_qPoints.at(qBin);
                    const double r = m_rPoints.at(rBin);
                    const double cosTheta = m_ctPoints.at(ctBin);
                    if (q < 2 || r < 1)
                    { 
                        m_grid(qBin,rBin,ctBin) = 0;
                    }
                    else
                    {
                        m_grid(qBin,rBin,ctBin) = m_waveFunction.GetPsiSquared(q,r,(cosTheta < -1) ? -1 : (cosTheta > 1) ? 1 : cosTheta); // limiting the value to stay within (-1,1)
                    }
                    
                }
    }

} // namespace JJCorrFitter
