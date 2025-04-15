/**
 * @file DoubleGaussian1D.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Defines the one-dimensional double-Gaussian source function
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

 #ifndef DoubleGaussian1D_hxx
 #define DoubleGaussian1D_hxx

 #include <vector>

 #include "SourceFunctionImpl.hxx"

 namespace JJCorrFitter
 {
     class DoubleGaussian1D : public SourceFunctionImpl
     {
         private:
             double m_invariantRadius1,m_invariantRadius2;

         public:
             DoubleGaussian1D(/* args */);
             ~DoubleGaussian1D() = default;
             DoubleGaussian1D(const DoubleGaussian1D&) = default;
             DoubleGaussian1D& operator=(const DoubleGaussian1D&) = default;
             DoubleGaussian1D(DoubleGaussian1D&&) noexcept = default;
             DoubleGaussian1D& operator=(DoubleGaussian1D&&) noexcept = default;

             void SetParameters(const std::vector<double> &pars);
             [[nodiscard]] double GetValue(float rStar) const noexcept;
             [[nodiscard]] double GetValue(float rOut, float rSide, float rLong) const noexcept {std::cerr << "DoubleGaussian1D::GetValue - This is a deleted function" << std::endl; return 0;}
     };

 } // namespace JJCorrFitter
 

#endif