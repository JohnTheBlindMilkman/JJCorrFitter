/**
 * @file SourceFunction3D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Defines the three-dimensional Gaussian source function
 * @version 0.1.0
 * @date 2024-10-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SourceFunction3D_hxx
    #define SourceFunction3D_hxx

    #include <vector>

    #include "SourceFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class SourceFunction3D : public SourceFunctionImpl
        {
            private:
                double m_outRadius, m_sideRadius, m_longRadius;

            public:
                SourceFunction3D(/* args */);
                ~SourceFunction3D() = default;
                SourceFunction3D(const SourceFunction3D&) = default;
                SourceFunction3D& operator=(const SourceFunction3D&) = default;
                SourceFunction3D(SourceFunction3D&&) noexcept = default;
                SourceFunction3D& operator=(SourceFunction3D&&) noexcept = default;

                void SetParameters(const std::vector<double> &pars);
                [[nodiscard]] double GetValue(float rStar) const noexcept {std::cerr << "SourceFunction3D::GetValue - This is a deleted function" << std::endl; return 0;}
                [[nodiscard]] double GetValue(float rOut, float rSide, float rLong) const noexcept;
        };

    } // namespace JJCorrFitter
    

#endif