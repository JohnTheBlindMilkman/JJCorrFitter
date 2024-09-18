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
            private:
                std::unique_ptr<TH1> m_dataToFit;
                std::unique_ptr<CorrelationFunctionImpl> m_corrFunction;
                std::unique_ptr<ROOT::Math::Minimizer> m_minimiser;
                std::unique_ptr<LikelihoodImpl> m_likelyhoodTest;

                void PrintInfo() const;

            public:
                Fitter() = delete;
                Fitter(std::unique_ptr<TH1> &&data, std::unique_ptr<CorrelationFunctionImpl> &&function, std::unique_ptr<ROOT::Math::Minimizer> &&minimiser ,std::unique_ptr<LikelihoodImpl> &&test);
                ~Fitter() = default;
                Fitter(const Fitter&) = delete;
                Fitter& operator=(const Fitter&) = delete;
                Fitter(Fitter&&) noexcept;
                Fitter& operator=(Fitter&&) noexcept;

                void SetHistogram(std::unique_ptr<TH1> &&data) noexcept;
                void SetCorrelationFunction(std::unique_ptr<CorrelationFunctionImpl> &&function) noexcept;
                void SetMinimiser(std::unique_ptr<ROOT::Math::Minimizer> &&minimiser) noexcept;
                void SetGoodnessTest(std::unique_ptr<LikelihoodImpl> &&test) noexcept;
                void SetParameter(int ival, const std::string &name, float start, float step, float min, float max);
                void SetParameter(int ival, const std::string &name, float start);
                bool Fit();
        };

    } // namespace JJCorrFitter

#endif