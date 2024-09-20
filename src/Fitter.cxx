#include "Fitter.hxx"

namespace JJCorrFitter
{
    Fitter::Fitter(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser ,std::unique_ptr<LikelihoodImpl> &&test) :
    m_currentParNumber(0), m_minimiser(std::move(minimiser)), m_likelyhoodTest(std::move(test)),
    m_parCounter{{ParType::Generic,{}},{ParType::Source,{}},{ParType::Interaction,{}}}
    {
        if (m_minimiser == nullptr)
            throw std::runtime_error("Fitter::Fitter - Error: cannot create minimizer");

        m_minimiser->SetMaxFunctionCalls(1000000);
        m_minimiser->SetTolerance(0.001);
        m_minimiser->SetPrintLevel(1);
    }

    Fitter::Fitter(Fitter &&other) noexcept : 
    m_minimiser(std::exchange(other.m_minimiser,nullptr)), m_likelyhoodTest(std::exchange(other.m_likelyhoodTest,nullptr)),
    m_parCounter(std::move(other.m_parCounter))
    {}

    Fitter& Fitter::operator=(Fitter &&other) noexcept
    {
        m_minimiser = std::exchange(other.m_minimiser,nullptr);
        m_likelyhoodTest = std::exchange(other.m_likelyhoodTest,nullptr);
        m_parCounter = std::move(other.m_parCounter);
    }

    void Fitter::SetMinimiser(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser) noexcept
    {
        m_minimiser = std::move(minimiser);
    }

    void Fitter::SetGoodnessTest(std::unique_ptr<LikelihoodImpl> &&test) noexcept
    {
        m_likelyhoodTest = std::move(test);
    }
    
    void Fitter::SetParameter(ParType type, const std::string &name, float start, float step, float min, float max)
    {
        //m_minimiser->SetLimitedVariable(m_currentParNumber,name,start,step,min,max);
        m_parCounter[type].push_back(m_currentParNumber);
        ++m_currentParNumber;
    }

    void Fitter::SetParameter(ParType type, const std::string &name, float start)
    {
        //m_minimiser->SetFixedVariable(m_currentParNumber,name,start);
        m_parCounter[type].push_back(m_currentParNumber);
        ++m_currentParNumber;
    }

    bool Fitter::Fit()
    {
        ROOT::Math::Functor f(m_likelyhoodTest->GetObjectiveFunction(m_parCounter.at(ParType::Generic),m_parCounter.at(ParType::Source),m_parCounter.at(ParType::Interaction)),m_likelyhoodTest->GetNParams());
        m_minimiser->SetFunction(f);

        m_minimiser->SetLimitedVariable(0,"N",50,1,10,100);
        m_minimiser->SetFixedVariable(1,"Lambda",1);
        m_minimiser->SetLimitedVariable(2,"Rinv",3,0.1,1,6); // I need to fix this somehow. Par limits can be set only after we set the function

        return m_minimiser->Minimize(); // TH1 is being overriden and produces waring and potentially memory leaks. Maybe we switch to vectors?
    }

    void Fitter::PrintInfo() const
    {
        std::cout << Config::projectName << "\n";
        std::cout << Config::projectVersion << "\n";
    }

} // namespace JJCorrFitter
