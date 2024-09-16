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

    namespace JJCorrFitter
    {
        class CorrelationFunctionImpl
        {
            private:

            public:
                CorrelationFunctionImpl(/* args */) = default;
                virtual ~CorrelationFunctionImpl(){}
                CorrelationFunctionImpl(const CorrelationFunctionImpl&) = default;
                CorrelationFunctionImpl& operator=(const CorrelationFunctionImpl&) = default;
                CorrelationFunctionImpl(CorrelationFunctionImpl&&) noexcept = default;
                CorrelationFunctionImpl& operator=(CorrelationFunctionImpl&&) noexcept = default;
                [[nodiscard]] virtual double Evaluate() = 0;
        };

    } // namespace JJCorrFitter
    

#endif