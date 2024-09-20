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
            public: 
                enum class ParType{Generic,Source,Interaction};

            private:
                std::size_t m_currentParNumber;
                std::unique_ptr<ROOT::Math::Minimizer> m_minimiser;
                std::unique_ptr<LikelihoodImpl> m_likelyhoodTest;
                std::unordered_map<ParType,std::vector<std::size_t> > m_parCounter;

                void PrintInfo() const;

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
                void SetParameter(ParType type, const std::string &name, float start, float step, float min, float max);
                void SetParameter(ParType type, const std::string &name, float start);
                bool Fit();
        };

    } // namespace JJCorrFitter

#endif