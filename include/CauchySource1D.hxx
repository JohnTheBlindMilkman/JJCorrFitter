/**
 * @file CauchySource1D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Defines the one-dimensional Cauchy source function
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

 #ifndef CauchySource1D_hxx
 #define CauchySource1D_hxx

 #include <vector>

 #include "SourceFunctionImpl.hxx"

 namespace JJCorrFitter
 {
     class CauchySource1D : public SourceFunctionImpl
     {
         private:
             double m_invariantRadius;

         public:
             CauchySource1D(/* args */);
             ~CauchySource1D() = default;
             CauchySource1D(const CauchySource1D&) = default;
             CauchySource1D& operator=(const CauchySource1D&) = default;
             CauchySource1D(CauchySource1D&&) noexcept = default;
             CauchySource1D& operator=(CauchySource1D&&) noexcept = default;

             void SetParameters(const std::vector<double> &pars);
             [[nodiscard]] double GetValue(float rStar) const noexcept;
             [[nodiscard]] double GetValue(float rOut, float rSide, float rLong) const noexcept {std::cerr << "CauchySource1D::GetValue - This is a deleted function" << std::endl; return 0;}
     };

 } // namespace JJCorrFitter
 

#endif