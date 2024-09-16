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
    #include "TH1D.h"

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
                float m_cosTheta, m_minKStar,m_maxKStar;
                int m_nPoints;
                const std::string m_histogramName{"hCF1D"}, m_histogramTitle{"One-dimensional correlation function"};
                std::vector<double> m_correlationFunctionPoints;

                [[nodiscard]] double CalculatePoint();
                [[nodiscard]] double CalculatePoint(float kStar);

            public:
                CorrelationFunction1D(/* args */) = delete;
                CorrelationFunction1D(float kStarMin, float kStarMax, int nPoints);
                ~CorrelationFunction1D() = default;
                CorrelationFunction1D(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D& operator=(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D(CorrelationFunction1D&&) noexcept = default;
                CorrelationFunction1D& operator=(CorrelationFunction1D&&) noexcept = default;

                [[nodiscard]] std::unique_ptr<TH1> Evaluate();
                void SetIntegrationRange(float rStarMin, float rStarMax) noexcept;
                void SetParameters() noexcept;
                void SetParameters(float rInv) noexcept;
        };
        
    } // namespace JJCorrFitter
    

#endif