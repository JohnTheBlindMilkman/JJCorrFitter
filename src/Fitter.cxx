#include "Fitter.hxx"

namespace JJCorrFitter
{
    Fitter::Fitter(std::unique_ptr<TH1> &&data, std::unique_ptr<CorrelationFunctionImpl> &&function, std::unique_ptr<ROOT::Math::Minimizer> &&minimiser ,std::unique_ptr<LikelihoodImpl> &&test) :
    m_dataToFit(std::move(data)), m_corrFunction(std::move(function)), m_minimiser(std::move(minimiser)), m_likelyhoodTest(std::move(test))
    {
        if (m_minimiser == nullptr)
            throw std::runtime_error("Fitter::Fitter - Error: cannot create minimizer");

        m_minimiser->SetMaxFunctionCalls(1000000);
        m_minimiser->SetTolerance(0.001);
        m_minimiser->SetPrintLevel(1);
    }

    Fitter::Fitter(Fitter &&other) noexcept : 
    m_dataToFit(std::exchange(other.m_dataToFit,nullptr)), m_corrFunction(std::exchange(other.m_corrFunction,nullptr)), 
    m_minimiser(std::exchange(other.m_minimiser,nullptr)), m_likelyhoodTest(std::exchange(other.m_likelyhoodTest,nullptr))
    {}

    Fitter& Fitter::operator=(Fitter &&other) noexcept
    {
        m_dataToFit = std::exchange(other.m_dataToFit,nullptr);
        m_corrFunction = std::exchange(other.m_corrFunction,nullptr);
        m_minimiser = std::exchange(other.m_minimiser,nullptr);
        m_likelyhoodTest = std::exchange(other.m_likelyhoodTest,nullptr);
    }

    void Fitter::SetHistogram(std::unique_ptr<TH1> &&data) noexcept
    {
        m_dataToFit = std::move(data);
    }

    void Fitter::SetCorrelationFunction(std::unique_ptr<CorrelationFunctionImpl> &&function) noexcept
    {
        m_corrFunction = std::move(function);
    }

    void Fitter::SetMinimiser(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser) noexcept
    {
        m_minimiser = std::move(minimiser);
    }

    void Fitter::SetGoodnessTest(std::unique_ptr<LikelihoodImpl> &&test) noexcept
    {
        m_likelyhoodTest = std::move(test);
    }
    
    void Fitter::SetParameter(int ival, const std::string &name, float start, float step, float min, float max)
    {
        m_minimiser->SetLimitedVariable(ival,name,start,step,min,max);
    }

    void Fitter::SetParameter(int ival, const std::string &name, float start)
    {
        m_minimiser->SetFixedVariable(ival,name,start);
    }

    bool Fitter::Fit()
    {
        ROOT::Math::Functor f(m_likelyhoodTest->GetTestFunction(),3);
        m_minimiser->SetFunction(f);
        return m_minimiser->Minimize();
    }

    void Fitter::PrintInfo() const
    {
        std::cout << Config::projectName << "\n";
        std::cout << Config::projectVersion << "\n";
    }

} // namespace JJCorrFitter
