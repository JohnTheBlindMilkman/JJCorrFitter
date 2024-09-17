/**
 * @file Fitter1D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Main class for performing 1D fitting of correlation functions
 * @version 0.1.0
 * @date 2024-09-17
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef Fitter1D_hxx
    #define Fitter1D_hxx

    #include <algorithm>

    #include "TH1.h"

    #include "CorrelationFunction1D.hxx"
    #include "FitterImpl.hxx"

    namespace JJCorrFitter
    {
        struct Fitter1DHelper
        {
            unsigned int number;
            std::string name;
            float step;
            bool isSet = false;
        };

        class Fitter1D : public FitterImpl
        {
            public:
                enum class Parameter{Norm, Lambda, Radius};

            private:
                std::map<Parameter,Fitter1DHelper> m_parameters;
                CorrelationFunction1D m_corrFunction;
                std::unique_ptr<ROOT::Math::Minimizer> m_minimiser;
                std::unique_ptr<TH1> m_dataToFit;

                bool AllParamsAreSet(const std::map<Parameter,Fitter1DHelper> &pars) noexcept;

            public:
                Fitter1D(/* args */) = delete;
                explicit Fitter1D(std::unique_ptr<TH1> &&data);
                ~Fitter1D() = default;
                Fitter1D(const Fitter1D&) = delete;
                Fitter1D& operator=(const Fitter1D&) = delete;
                Fitter1D(Fitter1D&&) noexcept;
                Fitter1D& operator=(Fitter1D&&) noexcept;

                bool Fit();
                void SetParameter(Parameter par, float start, float min, float max);
                void SetParameter(Parameter par, float start);
        };

    } // namespace JJCorrFitter
    

#endif