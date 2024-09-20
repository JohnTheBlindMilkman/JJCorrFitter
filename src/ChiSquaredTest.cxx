#include "ChiSquaredTest.hxx"

namespace JJCorrFitter
{
    ChiSquaredTest::ChiSquaredTest(std::unique_ptr<TH1> &&data, std::unique_ptr<CorrelationFunctionImpl> &&func)
    {
        m_dataToFit = std::move(data);
        m_corrFunc = std::move(func);
    }

    ChiSquaredTest::ChiSquaredTest(ChiSquaredTest &&other) noexcept
    {
        m_dataToFit = std::exchange(other.m_dataToFit,nullptr);
        m_corrFunc = std::exchange(other.m_corrFunc,nullptr);
    }

    ChiSquaredTest& ChiSquaredTest::operator=(ChiSquaredTest &&other) noexcept
    {
        m_dataToFit = std::exchange(other.m_dataToFit,nullptr);
        m_corrFunc = std::exchange(other.m_corrFunc,nullptr);

        return *this;
    }

    std::function<double (const double *)> ChiSquaredTest::GetObjectiveFunction(const std::vector<std::size_t> &corrFuncIndexes,const std::vector<std::size_t> &srcIndexes,const std::vector<std::size_t> &psiIndexes)
    {
        auto lambda = [=](const double *x)
        {
            m_corrFunc->SetParameters(this->SortParameters(x,corrFuncIndexes),this->SortParameters(x,srcIndexes),this->SortParameters(x,psiIndexes));
            return m_dataToFit->Chi2Test(m_corrFunc->Evaluate().release());
        };

        return lambda;
    }

} // namespace JJCorrFitter
