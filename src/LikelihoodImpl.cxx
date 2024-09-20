#include "LikelihoodImpl.hxx"

namespace JJCorrFitter
{
    void LikelihoodImpl::SetHistogram(std::unique_ptr<TH1> &&data) noexcept
    {
        m_dataToFit = std::move(data);
    }

    void LikelihoodImpl::SetCorrelationFunction(std::unique_ptr<CorrelationFunctionImpl> &&function) noexcept
    {
        m_corrFunc = std::move(function);
    }

    std::vector<double> LikelihoodImpl::SortParameters(const double *x, const std::vector<std::size_t> &indexArray) const
    {
        std::vector<double> tmp;
        for (const auto &index : indexArray)
            tmp.push_back(x[index]);

        return tmp;
    }
} // namespace JJCorrFitter
