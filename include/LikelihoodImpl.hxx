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
    #include <iostream>
    #include <chrono>
    #include <iomanip>

    #include "TH1.h"

    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class LikelihoodImpl
        {
            protected:
                static constexpr int m_rootHistogramFirstBin{1};

                std::unique_ptr<CorrelationFunctionImpl> m_corrFunc;
                std::unique_ptr<TH1> m_dataToFit;
                std::string_view m_likelihoodTestName,m_likelihoodResultName;

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
                [[nodiscard]] std::unique_ptr<TH1> GetHistogram() noexcept;
                [[nodiscard]] std::unique_ptr<TH1> GetCorrelationFunction() noexcept;
                [[nodiscard]] std::string_view GetLikelihoodTestName() const noexcept;
                [[nodiscard]] std::string_view GetLikelihoodResultName() const noexcept;
                [[nodiscard]] std::string_view GetCorrelationFunctionName() const noexcept;
                [[nodiscard]] std::string_view GetCorrelationSourceType() const noexcept;
                [[nodiscard]] std::string_view GetCorrelationInteractionTermType() const noexcept;
        };

            inline void LikelihoodImpl::SetHistogram(std::unique_ptr<TH1> &&data) noexcept {m_dataToFit = std::move(data);}
            inline void LikelihoodImpl::SetCorrelationFunction(std::unique_ptr<CorrelationFunctionImpl> &&function) noexcept {m_corrFunc = std::move(function);}
            inline std::unique_ptr<TH1> LikelihoodImpl::GetHistogram() noexcept {return std::move(m_dataToFit);}
            inline std::unique_ptr<TH1> LikelihoodImpl::GetCorrelationFunction() noexcept {return m_corrFunc->GetCorrelationFunction();}
            inline std::string_view LikelihoodImpl::GetLikelihoodTestName() const noexcept {return m_likelihoodTestName;}
            inline std::string_view LikelihoodImpl::GetLikelihoodResultName() const noexcept {return m_likelihoodResultName;}
            inline std::string_view LikelihoodImpl::GetCorrelationFunctionName() const noexcept {return m_corrFunc->GetFunctionName();}
            inline std::string_view LikelihoodImpl::GetCorrelationSourceType() const noexcept {return m_corrFunc->GetSourceType();}
            inline std::string_view LikelihoodImpl::GetCorrelationInteractionTermType() const noexcept {return m_corrFunc->GetInteractionTermType();}

    } // namespace JJCorrFitter
    
#endif