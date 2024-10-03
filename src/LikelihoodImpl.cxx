#include "LikelihoodImpl.hxx"

namespace JJCorrFitter
{
    std::vector<double> LikelihoodImpl::SortParameters(const double *x, const std::vector<std::size_t> &indexArray) const
    {
        std::vector<double> tmp;
        for (const auto &index : indexArray)
            tmp.push_back(x[index]);

        return tmp;
    }
} // namespace JJCorrFitter
