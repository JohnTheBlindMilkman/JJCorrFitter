/**
 * @file CorrelationFunction3D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief The correlation function object. Performs numerical integration of the Koonin-Pratt equation for 3D functions.
 * @version 0.1.0
 * @date 2024-09-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CorrelationFunction3D_hxx
    #define CorrelationFunction3D_hxx

    #include "Math/IntegratorMultiDim.h"
    #include "Math/Integrator.h"
    #include "TH3D.h"

    #include "gsl/gsl_sf_legendre.h"

    #include "CorrelationFunctionImpl.hxx"
    #include "Grid.hxx"

    namespace JJCorrFitter
    {
        class CorrelationFunction3D : public CorrelationFunctionImpl
        {
            private:
                double m_minKStar,m_maxKStar;
                int m_nPoints;
                std::string m_histogramName, m_histogramTitle;
                Grid<std::tuple<float,float,float> > m_kStarValues;
                Grid<double> m_correlationPoints, m_correlationErrors;

                [[nodiscard]] std::unique_ptr<TH3D> MakeHistogram(const Grid<double> &points, const Grid<double> &errors, const std::vector<double> &params);
                [[nodiscard]] Grid<std::tuple<float,float,float> > SetKStarPoints(float start, float stop, int nPoints);
                void NormaliseFunction([[maybe_unused]] std::vector<double> &points,[[maybe_unused]] std::vector<double> &errors) {}
                void NormaliseFunction(Grid<double> &points, Grid<double> &errors);
                [[nodiscard]] std::pair<double,double> CalculatePoint();
                [[nodiscard]] std::pair<double,double> CalculatePoint(float kOut,float kSide,float kLong);
                [[nodiscard]] constexpr double CalculateCosTheta(float kOut,float kSide,float kLong, float rOut,float rSide,float rLong) noexcept;
                [[nodiscard]] constexpr double CalculateModulus(float xOut,float xSide,float xLong) noexcept;

            public:
                CorrelationFunction3D(/* args */) = delete;
                CorrelationFunction3D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact);
                ~CorrelationFunction3D() = default;
                CorrelationFunction3D(const CorrelationFunction3D&) = delete;
                CorrelationFunction3D& operator=(const CorrelationFunction3D&) = delete;
                CorrelationFunction3D(CorrelationFunction3D&&) noexcept = default;
                CorrelationFunction3D& operator=(CorrelationFunction3D&&) noexcept = default;

                void SetBinning(const std::unique_ptr<TH1> &data, double minKstar = -1, double maxKstar = -1);
                void SetBinning(const std::string &name, const std::string &title, int nPoints, double minKStar, double maxKstar);
                [[nodiscard]] std::unique_ptr<TH1> Evaluate();
                [[nodiscard]] std::unique_ptr<TH1> EvaluateAtEdges();
                void SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars);
                [[nodiscard]] std::size_t GetNParams() const noexcept;
                [[nodiscard]] std::unique_ptr<TH1> GetCorrelationFunction();
        };

        [[nodiscard]] inline std::size_t CorrelationFunction3D::GetNParams() const noexcept {return m_numberOfParams;} 
        [[nodiscard]] inline std::unique_ptr<TH1> CorrelationFunction3D::GetCorrelationFunction() 
            {return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);}
        
    } // namespace JJCorrFitter
    

#endif