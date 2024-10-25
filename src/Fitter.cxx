#include "Fitter.hxx"

namespace JJCorrFitter
{
    Fitter::Fitter(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser ,std::unique_ptr<LikelihoodImpl> &&test) :
    m_currentParNumber(0), m_minimiser(std::move(minimiser)), m_likelyhoodTest(std::move(test)),
    m_parMap{{ParType::Generic,{}},{ParType::Source,{}},{ParType::Interaction,{}}}, m_valuesAtMinimum({}), m_errorsAtMinimum({}),
    m_tolerance(0.001), m_minFunctionValue(0.0), m_maxFunctionCalls(1000000),m_printLevel(0)
    {
        if (m_minimiser == nullptr)
            throw std::runtime_error("Fitter::Fitter - Error: cannot create minimizer");

        m_minimiser->SetMaxFunctionCalls(m_maxFunctionCalls);
        m_minimiser->SetTolerance(m_tolerance);
        m_minimiser->SetPrintLevel(m_printLevel);

        PrintHello();
    }

    Fitter::Fitter(Fitter &&other) noexcept : 
    m_minimiser(std::exchange(other.m_minimiser,nullptr)), m_likelyhoodTest(std::exchange(other.m_likelyhoodTest,nullptr)),
    m_parMap(std::move(other.m_parMap))
    {}

    Fitter& Fitter::operator=(Fitter &&other) noexcept
    {
        m_minimiser = std::exchange(other.m_minimiser,nullptr);
        m_likelyhoodTest = std::exchange(other.m_likelyhoodTest,nullptr);
        m_parMap = std::move(other.m_parMap);

        return *this;
    }
    
    void Fitter::SetParameter(ParType type, const std::string &name, float start, float step, float min, float max)
    {
        m_parMap[type].push_back({m_currentParNumber,name,start,step,min,max,false});
        ++m_currentParNumber;
    }

    void Fitter::SetParameter(ParType type, const std::string &name, float start)
    {
        m_parMap[type].push_back({m_currentParNumber,name,start,0,0,0,true});
        ++m_currentParNumber;
    }

    bool Fitter::Fit()
    {
        ROOT::Math::Functor f(
            m_likelyhoodTest->GetObjectiveFunction(
                GetIndices(m_parMap.at(ParType::Generic)),
                GetIndices(m_parMap.at(ParType::Source)),
                GetIndices(m_parMap.at(ParType::Interaction))),
            m_likelyhoodTest->GetNParams());

        m_minimiser->SetFunction(f);
        PassParametersToMinimiser(m_parMap);

        bool result = m_minimiser->Minimize();
        std::unique_ptr<const double> minVals(m_minimiser->X());
        m_valuesAtMinimum = PassPointerArrayToVector(m_minimiser->X(),m_likelyhoodTest->GetNParams());
        m_errorsAtMinimum = PassPointerArrayToVector(m_minimiser->Errors(),m_likelyhoodTest->GetNParams());

        PrintFitResult();

        return result;
    }

    void Fitter::PrintInfo() const
    {
        ROOT::Math::MinimizerOptions opts = m_minimiser->Options();
        std::cout << "-----===== Fitter Setup =====-----\n";
        std::cout << "Max function calls: " << m_maxFunctionCalls << "\n";
        std::cout << "Tolerance: " << m_tolerance << "\n";
        std::cout << "Minimizer print level: " << m_printLevel << "\n";
        std::cout << "Minimizer type: " << opts.MinimizerType() << "\n";
        std::cout << "Minimizer algorithm: " << opts.MinimizerAlgorithm() << "\n";
        std::cout << "Likelihood test: " << m_likelyhoodTest->GetLikelihoodTestName() << "\n";
        std::cout << "Correlation function type: " << m_likelyhoodTest->GetCorrelationFunctionName() << "\n";
        std::cout << "Source function: " << m_likelyhoodTest->GetCorrelationSourceType() << "\n";
        std::cout << "Squared wavefunction modulus: " << m_likelyhoodTest->GetCorrelationInteractionTermType() << " \n";
        std::cout << "-----========================-----\n";
    }

    void Fitter::PrintFitResult() const
    {
        std::cout << "-----===== Fit  Results =====-----\n";

        for (const auto &elem : m_parMap)
        {
            for (const auto &param : elem.second)
            {
                std::string_view str = param.isFixed ? "(fixed)" : " +/- " + std::to_string(m_errorsAtMinimum.at(param.number));
                std::cout << "Parameter " << param.name << ": " << m_valuesAtMinimum.at(param.number) << " " << str << "\n";
            }
        }
        std::cout << "----------------------------------\n";
        std::cout << m_likelyhoodTest->GetLikelihoodResultName() << ": " << m_minimiser->MinValue() << "\n";
        std::cout << "No. calls: " << m_minimiser->NCalls() << "\n";
        std::cout << "-----========================-----\n";
    }

    std::vector<std::size_t> Fitter::GetIndices(const std::vector<ParInfo> &vec) const noexcept
    {
        std::vector<std::size_t> tmp;
        tmp.reserve(vec.size());

        for (const auto &elem : vec)
        {
            tmp.push_back(elem.number);
        }

        return tmp;
    }

    std::vector<double> Fitter::PassPointerArrayToVector(const double *arr, std::size_t size)
    {
        std::vector<double> vec(size,0);
        if (arr != nullptr)
        {
            for (std::size_t i= 0; i < size; ++i)
                vec.at(i) = arr[i];
        }

        return vec;
    }

    void Fitter::PassParametersToMinimiser(const std::map<ParType,std::vector<ParInfo> > &params)
    {
        for (const auto &[key,vec] : params)
            for (const auto &elem : vec)
            {
                if (elem.isFixed)
                {
                    m_minimiser->SetFixedVariable(elem.number,elem.name,elem.start);
                }
                else
                {
                    m_minimiser->SetLimitedVariable(elem.number,elem.name,elem.start,elem.step,elem.min,elem.max);
                }
            }
    }

    void Fitter::PrintHello() const
    {
        std::cout << "-----===== "  << Config::m_projectName << " =====-----\n";
        std::cout << "Author: Jedrzej Kolas\n";
        std::cout << "Version: " << Config::m_projectVersion << "\n";
        std::cout << "Correlation function fitter for p-p.\n";
        std::cout << "Use at your own risk!\n";
        std::cout << "-----========================-----\n";
    }

} // namespace JJCorrFitter
