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

    #include <memory>

    #include "Math/Math.h"

    namespace JJCorrFitter
    {
        class SourceFunctionImpl
        {
            private:
                /* data */
            public:
                SourceFunctionImpl(/* args */) = default;
                virtual ~SourceFunctionImpl(){}
                SourceFunctionImpl(const SourceFunctionImpl&) = default;
                SourceFunctionImpl& operator=(const SourceFunctionImpl&) = default;
                SourceFunctionImpl(SourceFunctionImpl&&) noexcept = default;
                SourceFunctionImpl& operator=(SourceFunctionImpl&&) noexcept = default;

                virtual void SetParameters(float rInv) noexcept = 0;
                [[nodiscard]] virtual double GetValue(float rStar) const noexcept = 0;
        };
        
    } // namespace JJCorrFitter
    

#endif