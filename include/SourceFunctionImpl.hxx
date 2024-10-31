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
    #include <iostream>

    #include "Math/Math.h"

    namespace JJCorrFitter
    {
        class SourceFunctionImpl
        {
            private:
                /* data */

            protected:
                std::size_t m_numberOfParams;
                std::string_view m_sourceFunctionName;

            public:
                SourceFunctionImpl(/* args */) = default;
                virtual ~SourceFunctionImpl(){}
                SourceFunctionImpl(const SourceFunctionImpl&) = default;
                SourceFunctionImpl& operator=(const SourceFunctionImpl&) = default;
                SourceFunctionImpl(SourceFunctionImpl&&) noexcept = default;
                SourceFunctionImpl& operator=(SourceFunctionImpl&&) noexcept = default;

                virtual void SetParameters(const std::vector<double> &pars) = 0;
                [[nodiscard]] virtual double GetValue(float rStar) const noexcept = 0;
                [[nodiscard]] virtual double GetValue(float rOut, float rSide, float rLong) const noexcept = 0;
                [[nodiscard]] std::size_t GetNParams() const noexcept;
                [[nodiscard]] std::string_view GetSourceFunctionName() const noexcept;
        };

        inline std::size_t SourceFunctionImpl::GetNParams() const noexcept {return m_numberOfParams;}
        inline std::string_view SourceFunctionImpl::GetSourceFunctionName() const noexcept {return m_sourceFunctionName;}
        
    } // namespace JJCorrFitter
    

#endif