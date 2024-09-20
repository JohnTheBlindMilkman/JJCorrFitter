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
    #include <functional>

    #include "TH1.h"

    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class LikelihoodImpl
        {
            protected:
                std::unique_ptr<CorrelationFunctionImpl> m_corrFunc;
                std::unique_ptr<TH1> m_dataToFit;

                std::vector<double> SortParameters(const double *x, const std::vector<std::size_t> &indexArray) const;

            public:
                LikelihoodImpl(/* args */) = default;
                virtual ~LikelihoodImpl() = default;
                LikelihoodImpl(const LikelihoodImpl&) = delete;
                LikelihoodImpl& operator=(const LikelihoodImpl&) = delete;
                LikelihoodImpl(LikelihoodImpl&&) noexcept = default;
                LikelihoodImpl& operator=(LikelihoodImpl&&) noexcept = default;

                void SetHistogram(std::unique_ptr<TH1> &&data) noexcept;
                void SetCorrelationFunction(std::unique_ptr<CorrelationFunctionImpl> &&function) noexcept;
                [[nodiscard]] virtual std::function<double (const double *)> GetObjectiveFunction(const std::vector<std::size_t> &corrFuncIndexes,const std::vector<std::size_t> &srcIndexes,const std::vector<std::size_t> &psiIndexes) = 0;
                [[nodiscard]] virtual std::size_t GetNParams() const = 0;
        };
        
    } // namespace JJCorrFitter
    
#endif