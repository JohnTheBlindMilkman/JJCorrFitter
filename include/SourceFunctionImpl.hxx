/**
 * @file SourceFunctionImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Defines the Gaussian source function for a given dimension
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SourceFunctionImpl_hxx
    #define SourceFunctionImpl_hxx

    #include "Math/Math.h"

    namespace JJCorrFitter
    {
        class SourceFunctionImpl
        {
            private:
                /* data */
            public:
                SourceFunctionImpl(/* args */);
                virtual ~SourceFunctionImpl(){}
                SourceFunctionImpl(const SourceFunctionImpl&) = delete;
                SourceFunctionImpl& operator=(const SourceFunctionImpl&) = delete;
                SourceFunctionImpl(SourceFunctionImpl&&) noexcept = default;
                SourceFunctionImpl& operator=(SourceFunctionImpl&&) noexcept = default;
                [[nodiscard]] virtual double GetValue(float rStar, float rInv) const noexcept = 0;
                [[nodiscard]] virtual double GetValue(float rOut, float rSide, float rLong, float Rout, float Rside, float Rlong) const noexcept = 0;
        };
        
    } // namespace JJCorrFitter
    

#endif