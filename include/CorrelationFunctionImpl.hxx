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

    namespace JJCorrFitter
    {
        class CorrelationFunctionImpl
        {
            private:

            protected:
                float m_minRStar, m_MaxRStar;
                [[nodiscard]] virtual double CalculatePoint() = 0;

            public:
                CorrelationFunctionImpl(/* args */) : m_minRStar(0.f), m_MaxRStar(20.f) {}
                virtual ~CorrelationFunctionImpl(){}
                CorrelationFunctionImpl(const CorrelationFunctionImpl&) = default;
                CorrelationFunctionImpl& operator=(const CorrelationFunctionImpl&) = default;
                CorrelationFunctionImpl(CorrelationFunctionImpl&&) noexcept = default;
                CorrelationFunctionImpl& operator=(CorrelationFunctionImpl&&) noexcept = default;

                [[nodiscard]] virtual std::unique_ptr<TH1> Evaluate() = 0;
                virtual void SetIntegrationRange(float rStarMin, float rStarMax) noexcept = 0;
                virtual void SetParameters() noexcept = 0;
                [[nodiscard]] virtual int GetNParams() const = 0;
        };

    } // namespace JJCorrFitter
    
#endif