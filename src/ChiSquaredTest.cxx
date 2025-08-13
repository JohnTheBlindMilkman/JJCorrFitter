#include "ChiSquaredTest.hxx"

namespace JJCorrFitter
{
    ChiSquaredTest::ChiSquaredTest(std::unique_ptr<TH1> &&data, std::unique_ptr<CorrelationFunctionImpl> &&func)
    {
        m_likelihoodTestName = "Chi-squared";
        m_likelihoodResultName = "Chi2/ndf";
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
            m_corrFunc->SetParameters(SortParameters(x,corrFuncIndexes),SortParameters(x,srcIndexes),SortParameters(x,psiIndexes));
            double chi2Ndf = CalculateChi2(m_dataToFit,m_corrFunc->Evaluate());

            return chi2Ndf;
        };

        return lambda;
    }

    double ChiSquaredTest::CalculateChi2(const std::unique_ptr<TH1> &data, const std::unique_ptr<TH1> &model)
    {
        // if (std::abs(data->GetBinWidth(m_rootHistogramFirstBin) - model->GetBinWidth(m_rootHistogramFirstBin)) > std::sqrt(std::numeric_limits<double>::epsilon()))
        // {
        //     throw std::logic_error(
        //         "ChiSquaredTest::CalculateChi2 - Bin width differs between data and models. \nData bin width is " + 
        //         std::to_string(data->GetBinWidth(m_rootHistogramFirstBin)) + 
        //         "\nModel bin width is " + 
        //         std::to_string(model->GetBinWidth(m_rootHistogramFirstBin))
        //     );
        // }
        
        const int nBins = model->GetNbinsX();
        double result = 0;
        int ndf = 0;
        for (int bin = 1; bin <= nBins; ++bin)
        {
            double modelVal = model->GetBinContent(bin);
            double dataVal = data->GetBinContent(data->FindBin(model->GetBinCenter(bin)));
            double dataErr = data->GetBinError(bin);

            if (dataVal > 0 && modelVal > 0 && dataErr > 0)
            {
                result += (dataVal - modelVal) * (dataVal - modelVal) / (dataErr * dataErr);
                ++ndf;
            }
        }

        return result / (ndf - 1);
    }

} // namespace JJCorrFitter
