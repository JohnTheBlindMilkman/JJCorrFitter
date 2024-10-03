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

    #include "Math/IntegratorMultiDim.h"
    #include "Math/Integrator.h"
    #include "TH1D.h"

    #include "gsl/gsl_sf_legendre.h"

    #include "boost/math/quadrature/gauss.hpp"
    #include "boost/math/quadrature/gauss_kronrod.hpp"

    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class CorrelationFunction1D : public CorrelationFunctionImpl
        {
            private:
                float m_minKStar,m_maxKStar;
                int m_nPoints;
                std::string m_histogramName, m_histogramTitle;
                std::vector<float> m_kStarValues;

                [[nodiscard]] std::unique_ptr<TH1D> MakeHistogram(const std::vector<double> &points, const std::vector<double> &errors, const std::vector<double> &params);
                [[nodiscard]] std::vector<float> SetKStarPoints(float start, float stop, int nPoints);
                void NormaliseFunction(std::vector<double> &points, std::vector<double> &errors);
                [[nodiscard]] std::pair<double,double> CalculatePoint();
                [[nodiscard]] std::pair<double,double> CalculatePoint(float kStar);

            public:
                CorrelationFunction1D(/* args */) = delete;
                CorrelationFunction1D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact);
                ~CorrelationFunction1D() = default;
                CorrelationFunction1D(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D& operator=(const CorrelationFunction1D&) = delete;
                CorrelationFunction1D(CorrelationFunction1D&&) noexcept = default;
                CorrelationFunction1D& operator=(CorrelationFunction1D&&) noexcept = default;

                void SetBinning(const std::unique_ptr<TH1> &data, float minKstar = -1, float maxKstar = -1);
                void SetBinning(const std::string &name, const std::string &title, int nPoints, float minKStar, float maxKstar);
                [[nodiscard]] std::unique_ptr<TH1> Evaluate();
                void SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars);
                [[nodiscard]] std::size_t GetNParams() const noexcept;
                [[nodiscard]] std::unique_ptr<TH1> GetCorrelationFunction();
        };

        [[nodiscard]] inline std::size_t CorrelationFunction1D::GetNParams() const noexcept {return m_numberOfParams;} 
        [[nodiscard]] inline std::unique_ptr<TH1> CorrelationFunction1D::GetCorrelationFunction() 
            {return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);}
        
    } // namespace JJCorrFitter
    

#endif