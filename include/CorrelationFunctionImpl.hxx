/**
 * @file CorrelationFunctionImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief The correlation function object. Performs numerical integration of the Koonin-Pratt equation.
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CorrelationFunctionImpl_hxx
    #define CorrelationFunctionImpl_hxx

    #include <memory>
    
    #include <TH1.h>

    #include "SourceFunctionImpl.hxx"
    #include "InteractionTermImpl.hxx"

    namespace JJCorrFitter
    {
        class CorrelationFunctionImpl
        {
            private:

            protected:
                std::size_t m_numberOfParams, m_totalNumberOfParams;
                std::unique_ptr<SourceFunctionImpl> m_sourceFunction;
                std::unique_ptr<InteractionTermImpl> m_interactionTerm;
                std::vector<double> m_corrFuncParams, m_sourceFunctionParams, m_interactionTermParams;

                [[nodiscard]] virtual double CalculatePoint() = 0;

            public:
                CorrelationFunctionImpl(/* args */) = default;
                virtual ~CorrelationFunctionImpl() = default;
                CorrelationFunctionImpl(const CorrelationFunctionImpl&) = delete;
                CorrelationFunctionImpl& operator=(const CorrelationFunctionImpl&) = delete;
                CorrelationFunctionImpl(CorrelationFunctionImpl&&) noexcept = default;
                CorrelationFunctionImpl& operator=(CorrelationFunctionImpl&&) noexcept = default;

                [[nodiscard]] virtual std::unique_ptr<TH1> Evaluate() = 0;
                virtual void SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars) = 0;
                [[nodiscard]] std::size_t GetNParams() const noexcept;
        };

        inline std::size_t CorrelationFunctionImpl::GetNParams() const noexcept {return m_totalNumberOfParams;}

    } // namespace JJCorrFitter
    
#endif