/**
 * @file LikelihoodImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Base class for impelmenting the goodness-of-fit function for the minimizer
 * @version 0.2.0
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef Likelihood_hxx
    #define Likelihood_hxx

    #include <memory>

    #include "TH1.h"

    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class LikelihoodImpl
        {
            protected:
                std::unique_ptr<CorrelationFunctionImpl> m_corrFunc;
                std::unique_ptr<TH1> m_dataToFit;

            public:
                LikelihoodImpl(/* args */) = default;
                virtual ~LikelihoodImpl() = default;
                LikelihoodImpl(const LikelihoodImpl&) = delete;
                LikelihoodImpl& operator=(const LikelihoodImpl&) = delete;
                LikelihoodImpl(LikelihoodImpl&&) noexcept = default;
                LikelihoodImpl& operator=(LikelihoodImpl&&) noexcept = default;

                [[nodiscard]] virtual std::function<double (const double *)> GetTestFunction() = 0;
                [[nodiscard]] virtual int GetNParams() const = 0;
        };
        
    } // namespace JJCorrFitter
    
#endif