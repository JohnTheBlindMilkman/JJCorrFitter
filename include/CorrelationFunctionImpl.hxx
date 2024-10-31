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
                static constexpr std::array<double,6> m_normalisationPoints{200,250,300,350,400,450};
                
                std::size_t m_numberOfParams, m_totalNumberOfParams;
                std::unique_ptr<SourceFunctionImpl> m_sourceFunction;
                std::unique_ptr<InteractionTermImpl> m_interactionTerm;
                std::vector<double> m_corrFuncParams, m_sourceFunctionParams, m_interactionTermParams;
                std::string_view m_correlationFunctionName;

                virtual void NormaliseFunction(std::vector<double> &points, std::vector<double> &errors) = 0;
                [[nodiscard]] virtual std::pair<double,double> CalculatePoint() = 0;

            public:
                CorrelationFunctionImpl(/* args */) = default;
                virtual ~CorrelationFunctionImpl() = default;
                CorrelationFunctionImpl(const CorrelationFunctionImpl&) = delete;
                CorrelationFunctionImpl& operator=(const CorrelationFunctionImpl&) = delete;
                CorrelationFunctionImpl(CorrelationFunctionImpl&&) noexcept = default;
                CorrelationFunctionImpl& operator=(CorrelationFunctionImpl&&) noexcept = default;

                virtual void SetBinning(const std::unique_ptr<TH1> &data, float minKstar = -1, float maxKstar = -1) = 0;
                virtual void SetBinning(const std::string &name, const std::string &title, int nPoints, float minKStar, float maxKstar) = 0;
                [[nodiscard]] virtual std::unique_ptr<TH1> Evaluate() = 0;
                virtual void SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars) = 0;
                [[nodiscard]] std::size_t GetNParams() const noexcept;
                [[nodiscard]] virtual std::unique_ptr<TH1> GetCorrelationFunction() = 0;
                [[nodiscard]] std::string_view GetFunctionName() const noexcept;
                [[nodiscard]] std::string_view GetSourceType() const noexcept;
                [[nodiscard]] std::string_view GetInteractionTermType() const noexcept;
        };

        inline std::size_t CorrelationFunctionImpl::GetNParams() const noexcept {return m_totalNumberOfParams;}
        inline std::string_view CorrelationFunctionImpl::GetFunctionName() const noexcept {return m_correlationFunctionName;}
        inline std::string_view CorrelationFunctionImpl::GetSourceType() const noexcept {return m_sourceFunction->GetSourceFunctionName();}
        inline std::string_view CorrelationFunctionImpl::GetInteractionTermType() const noexcept {return m_interactionTerm->GetInteractiontermName();}

    } // namespace JJCorrFitter
    
#endif