/**
 * @file FitterImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Main virtual class defining the fitter object
 * @version 0.1.0
 * @date 2024-09-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FitterImpl_hxx
    #define FitterImpl_hxx

    #include "Math/Minimizer.h"
    #include "Math/Factory.h"

    namespace JJCorrFitter
    {
        class FitterImpl
        {
            private:

            public:
                FitterImpl(){}
                virtual ~FitterImpl(){}
                FitterImpl(const FitterImpl&) = delete;
                FitterImpl& operator=(const FitterImpl&) = delete;
                FitterImpl(FitterImpl&&) noexcept = default;
                FitterImpl& operator=(FitterImpl&&) noexcept = default;

                virtual bool Fit() = 0;
        };

    } // namespace JJCorrFitter

#endif