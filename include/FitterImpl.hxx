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

    #include "CorrelationFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class FitterImpl
        {
            private:
                std::unique_ptr<CorrelationFunctionImpl> m_corrFunction;

            public:
                FitterImpl();
                virtual ~FitterImpl(){}
                virtual void Fit() = 0;
        };

    } // namespace JJCorrFitter

#endif