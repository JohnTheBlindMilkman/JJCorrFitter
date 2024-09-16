/**
 * @file CorrelationFunction1D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief The correlation function object. Performs numerical integration of the Koonin-Pratt equation for 1D functions.
 * @version 0.1.0
 * @date 2024-09-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CorrelationFunction1D_hxx
    #define CorrelationFunction1D_hxx

    #include "Math/Integrator.h"

    #include "CorrelationFunctionImpl.hxx"
    #include "InteractionTermSchrodinger.hxx"
    #include "SourceFunction1D.hxx"

    namespace JJCorrFitter
    {
        class CorrelationFunction1D : public CorrelationFunctionImpl
        {
            private:
                InteractionTermSchrodinger m_ppSchroed;
                SourceFunction1D m_source1D;

                static double IntegralExpression(double x);

            public:
                CorrelationFunction1D(/* args */);
                ~CorrelationFunction1D();
                CorrelationFunction1D(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D& operator=(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D(CorrelationFunction1D&&) noexcept = default;
                CorrelationFunction1D& operator=(CorrelationFunction1D&&) noexcept = default;
                [[nodiscard]] double Evaluate(int kStar, double rStar, double rInv);
        };
        
    } // namespace JJCorrFitter
    

#endif