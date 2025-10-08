#include "LogLikelihoodTest.hxx"

namespace JJCorrFitter
{
    LogLikelihoodTest::LogLikelihoodTest(std::unique_ptr<TH1> &&data, std::unique_ptr<TH1> &&signal, std::unique_ptr<TH1> &&background, std::unique_ptr<CorrelationFunctionImpl> &&model)
    {
        m_likelihoodTestName = "Log-likelihood";
        m_likelihoodResultName = "max log";
        m_dataToFit = std::move(data);
        m_dataSignal = std::move(signal);
        m_dataBackground = std::move(background);
        m_corrFunc = std::move(model);
    }

    LogLikelihoodTest::LogLikelihoodTest(LogLikelihoodTest &&other) noexcept
    {
        m_dataToFit = std::exchange(other.m_dataToFit,nullptr);
        m_dataSignal = std::exchange(other.m_dataSignal,nullptr);
        m_dataBackground = std::exchange(other.m_dataBackground,nullptr);
        m_corrFunc = std::exchange(other.m_corrFunc,nullptr);
    }

    LogLikelihoodTest& LogLikelihoodTest::operator=(LogLikelihoodTest &&other) noexcept
    {
        m_dataToFit = std::exchange(other.m_dataToFit,nullptr);
        m_dataSignal = std::exchange(other.m_dataSignal,nullptr);
        m_dataBackground = std::exchange(other.m_dataBackground,nullptr);
        m_corrFunc = std::exchange(other.m_corrFunc,nullptr);

        return *this;
    }

    std::function<double (const double *)> LogLikelihoodTest::GetObjectiveFunction(const std::vector<std::size_t> &corrFuncIndexes,const std::vector<std::size_t> &srcIndexes,const std::vector<std::size_t> &psiIndexes)
    {
        auto lambda = [this,corrFuncIndexes,srcIndexes,psiIndexes](const double *x)
        {
            m_corrFunc->SetParameters(SortParameters(x,corrFuncIndexes),SortParameters(x,srcIndexes),SortParameters(x,psiIndexes));
            return CalculateLogLikelihood(m_dataSignal,m_dataBackground,m_corrFunc->Evaluate());
        };

        return lambda;
    }

    double LogLikelihoodTest::CalculateLogLikelihood(const std::unique_ptr<TH1> &signal, const std::unique_ptr<TH1> &background, const std::unique_ptr<TH1> &model)
    {
        if (std::abs(signal->GetBinWidth(m_rootHistogramFirstBin) - model->GetBinWidth(m_rootHistogramFirstBin)) > std::sqrt(std::numeric_limits<Double_t>::epsilon()))
        {
            throw std::logic_error(
                "ChiSquaredTest::CalculateChi2 - Bin width differs between data and models. \nData bin width is " + 
                std::to_string(signal->GetBinWidth(m_rootHistogramFirstBin)) + 
                "\nModel bin width is " + 
                std::to_string(model->GetBinWidth(m_rootHistogramFirstBin))
            );
        }

        signal->Scale(1./signal->Integral("w"));
        background->Scale(1./background->Integral("w"));

        const int nBins = model->GetNbinsX();
        double result = 0, logA = 0, logB = 0, A = 0, B = 0, C = 0;
        for (int bin = 1; bin <= nBins; ++bin)
        {
            A = signal->GetBinContent(bin);
            B = background->GetBinContent(bin);
            C = model->GetBinContent(bin);
            // C = A / B; // maybe the normalisation is the problem?
            if (A > 0 && B > 0)
            {
                logA = (C * (A + B))/(A * (C + 1));
                logB = (A + B)/(B * (C + 1));
                // result += 2 * ((B * C) - (A * std::log(B * C))); // from Chris Ansen's code
                result += -2 * (A * std::log(logA) + B * std::log(logB)); // from Phys. Rev. C 66 (2002) 054906
            }
        }

        return result;
    }

} // namespace JJCorrFitter
