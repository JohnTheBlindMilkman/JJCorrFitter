/**
 * @file SourceFunction1D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Defines the one-dimensional Gaussian source function
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SourceFunction1D_hxx
    #define SourceFunction1D_hxx

    #include <vector>

    #include "SourceFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class SourceFunction1D : public SourceFunctionImpl
        {
            private:
                double m_invariantRadius;

            public:
                SourceFunction1D(/* args */);
                ~SourceFunction1D() = default;
                SourceFunction1D(const SourceFunction1D&) = default;
                SourceFunction1D& operator=(const SourceFunction1D&) = default;
                SourceFunction1D(SourceFunction1D&&) noexcept = default;
                SourceFunction1D& operator=(SourceFunction1D&&) noexcept = default;

                void SetParameters(const std::vector<double> &pars);
                [[nodiscard]] double GetValue(float rStar) const noexcept;
        };

    } // namespace JJCorrFitter
    

#endif