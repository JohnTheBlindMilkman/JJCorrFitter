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
        auto lambda = [this,corrFuncIndexes,srcIndexes,psiIndexes](const double *x)
        {
            m_corrFunc->SetParameters(SortParameters(x,corrFuncIndexes),SortParameters(x,srcIndexes),SortParameters(x,psiIndexes));
            double chi2Ndf = CalculateChi2(m_dataToFit,m_corrFunc->Evaluate());

            return chi2Ndf;
        };

        return lambda;
    }

    double ChiSquaredTest::CalculateChi2(const std::unique_ptr<TH1> &data, const std::unique_ptr<TH1> &model) const
    {
        static std::size_t count = 0;
        const auto dataWidth = data->GetXaxis()->GetBinWidth(m_rootHistogramFirstBin);
        const auto modelWidth = model->GetXaxis()->GetBinWidth(m_rootHistogramFirstBin);
        if (std::abs(dataWidth - modelWidth) > std::sqrt(std::max(dataWidth,modelWidth) * std::numeric_limits<double>::epsilon()))
        {
            throw std::runtime_error(
                "ChiSquaredTest::CalculateChi2 - Bin width differs between data and models. \nData bin width is " + 
                std::to_string(dataWidth) + 
                "\nModel bin width is " + 
                std::to_string(modelWidth)
            );
        }

        std::cout << "Count : " << ++count << "\n";

        std::string histClass(data->ClassName());
        if (histClass.find("TH1") != std::string::npos)
        {
            return CalculateChi2TH1(data,model);
        }
        else if (histClass.find("TH3") != std::string::npos)
        {
            return CalculateChi2TH3(data,model);
        }
        else
        {
            throw std::runtime_error("ChiSquaredTest::CalculateChi2 - given data histogram has unsuported type");
            return 0;
        }
    }

    double ChiSquaredTest::CalculateChi2TH1(const std::unique_ptr<TH1> &data, const std::unique_ptr<TH1> &model) const
    {
        const int nBins = model->GetNbinsX();
        double result = 0;
        int ndf = 0;
        for (int bin = 1; bin <= nBins; ++bin)
        {
            const double modelVal = model->GetBinContent(bin);
            const double dataVal = data->GetBinContent(data->FindBin(model->GetBinCenter(bin)));
            const double dataErr = data->GetBinError(bin);

            if (dataVal > 0 && modelVal > 0 && dataErr > 0)
            {
                result += (dataVal - modelVal) * (dataVal - modelVal) / (dataErr * dataErr);
                ++ndf;
            }
        }

        return result / (ndf - 1);
    }

    double ChiSquaredTest::CalculateChi2TH3(const std::unique_ptr<TH1> &data, const std::unique_ptr<TH1> &model) const
    {
        const int nBinsX = model->GetNbinsX();
        const int nBinsY = model->GetNbinsY();
        const int nBinsZ = model->GetNbinsZ();
        double result = 0;
        int ndf = 0;
        for (int binX = 1; binX <= nBinsX; ++binX)
            for (int binY = 1; binY <= nBinsY; ++binY)
                for (int binZ = 1; binZ <= nBinsZ; ++binZ)
                {
                    const double modelVal = model->GetBinContent(binX,binY,binZ);
                    const double dataVal = data->GetBinContent(
                        data->GetXaxis()->FindBin(model->GetXaxis()->GetBinCenter(binX)),
                        data->GetYaxis()->FindBin(model->GetYaxis()->GetBinCenter(binY)),
                        data->GetZaxis()->FindBin(model->GetZaxis()->GetBinCenter(binZ))
                    );
                    const double dataErr = data->GetBinError(binX,binY,binZ);

                    if (dataVal > 0 && modelVal > 0 && dataErr > 0)
                    {
                        result += (dataVal - modelVal) * (dataVal - modelVal) / (dataErr * dataErr);
                        ++ndf;
                    }
                }

        return result / (ndf - 1);
    }

} // namespace JJCorrFitter
