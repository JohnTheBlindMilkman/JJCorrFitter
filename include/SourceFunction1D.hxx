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

    #include "SourceFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class SourceFunction1D : public SourceFunctionImpl
        {
            private:
                /* data */
            public:
                SourceFunction1D(/* args */) = default;
                ~SourceFunction1D() = default;
                SourceFunction1D(const SourceFunction1D&) = delete;
                SourceFunction1D& operator=(const SourceFunction1D&) = delete;
                SourceFunction1D(SourceFunction1D&&) noexcept = default;
                SourceFunction1D& operator=(SourceFunction1D&&) noexcept = default;
                [[nodiscard]] double GetValue(float rStar, float rInv) const noexcept;
                [[nodiscard]] double GetValue([[maybe_unused]] float rOut, [[maybe_unused]] float rSide, [[maybe_unused]] float rLong, [[maybe_unused]] float Rout, [[maybe_unused]] float Rside, [[maybe_unused]] float Rlong) const noexcept {return 0;}
        };

    } // namespace JJCorrFitter
    

#endif