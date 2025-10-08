/**
 * @file Fitter.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Main class defining the fitter object
 * @version 0.1.0
 * @date 2024-09-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef Fitter_hxx
    #define Fitter_hxx

    #include <map>

    #include "Math/Minimizer.h"
    #include "Math/Factory.h"
    #include "Math/Functor.h"

    #include "Config.hxx"
    #include "LikelihoodImpl.hxx"
    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class Fitter
        {
            public: 
                enum class ParType{Generic,Source,Interaction};

            private:
                struct ParInfo
                {
                    std::size_t number;
                    std::string name;
                    float start,step,min,max;
                    bool isFixed;
                };

                std::size_t m_currentParNumber;
                std::unique_ptr<ROOT::Math::Minimizer> m_minimiser;
                std::unique_ptr<LikelihoodImpl> m_likelyhoodTest;
                std::map<ParType,std::vector<ParInfo> > m_parMap;
                std::vector<double> m_valuesAtMinimum, m_errorsAtMinimum;
                double m_tolerance, m_minFunctionValue;
                unsigned int m_maxFunctionCalls;
                int m_printLevel;

                [[nodiscard]] std::vector<std::size_t> GetIndices(const std::vector<ParInfo> &vec) const noexcept;
                [[nodiscard]] std::vector<double> PassPointerArrayToVector(const double *arr, std::size_t size);
                void PassParametersToMinimiser(const std::map<ParType,std::vector<ParInfo> > &params);
                void PrintHello() const;
                void PrintFitResult() const;

            public:
                Fitter() = delete;
                Fitter(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser ,std::unique_ptr<LikelihoodImpl> &&test);
                ~Fitter() = default;
                Fitter(const Fitter&) = delete;
                Fitter& operator=(const Fitter&) = delete;
                Fitter(Fitter&&) noexcept;
                Fitter& operator=(Fitter&&) noexcept;

                void SetMinimiser(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser) noexcept;
                void SetGoodnessTest(std::unique_ptr<LikelihoodImpl> &&test) noexcept;
                void SetMaxFunctionCalls(unsigned int maxCalls);
                void SetMaxIterations(unsigned maxIters) noexcept;
                void SetTolerance(double tol);
                void SetPrintLevel(int lvl);
                void SetParameter(ParType type, const std::string &name, float start, float step, float min, float max);
                void SetParameter(ParType type, const std::string &name, float start);
                void SetDataHistogram(std::unique_ptr<TH1> &&data) noexcept {m_likelyhoodTest->SetHistogram(std::move(data));}
                bool Fit();
                void PrintInfo() const;
                [[nodiscard]] std::unique_ptr<TH1> GetFitFunction() noexcept;
                [[nodiscard]] std::unique_ptr<TH1> GetDataHistogram() noexcept;
                [[nodiscard]] std::vector<double> GetFitParameterValues() const noexcept;
                [[nodiscard]] std::vector<double> GetFitParameterErrors() const noexcept;
                [[nodiscard]] int GetStatus() const noexcept {return m_minimiser->Status();}
                [[nodiscard]] double GetMinValue() const noexcept {return m_minimiser->MinValue();}
        };

        inline void Fitter::SetMinimiser(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser) noexcept {m_minimiser = std::move(minimiser);}
        inline void Fitter::SetGoodnessTest(std::unique_ptr<LikelihoodImpl> &&test) noexcept {m_likelyhoodTest = std::move(test);}
        inline void Fitter::SetMaxFunctionCalls(unsigned int maxCalls) {m_maxFunctionCalls = maxCalls; m_minimiser->SetMaxFunctionCalls(m_maxFunctionCalls);}
        inline void Fitter::SetMaxIterations(unsigned maxIters) noexcept {m_minimiser->SetMaxIterations(maxIters);}
        inline void Fitter::SetTolerance(double tol) {m_tolerance = tol; m_minimiser->SetTolerance(m_tolerance);}
        inline void Fitter::SetPrintLevel(int lvl) {m_printLevel = lvl; m_minimiser->SetPrintLevel(m_printLevel);}
        inline std::unique_ptr<TH1> Fitter::GetFitFunction() noexcept {return m_likelyhoodTest->GetCorrelationFunction();}
        inline std::unique_ptr<TH1> Fitter::GetDataHistogram() noexcept {return m_likelyhoodTest->GetHistogram();}
        inline std::vector<double> Fitter::GetFitParameterValues() const noexcept {return m_valuesAtMinimum;}
        inline std::vector<double> Fitter::GetFitParameterErrors() const noexcept {return m_errorsAtMinimum;}

    } // namespace JJCorrFitter

#endif