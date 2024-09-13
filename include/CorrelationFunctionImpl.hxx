/**
 * @file CorrelationFunctionImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief The correlation function object. Performs numerical integration of the Koonin-Pratt equation.
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CorrelationFunctionImpl_hxx
    #define CorrelationFunctionImpl_hxx

    #include <memory>

    #include "InteractionTermImpl.hxx"
    #include "SourceFunctionImpl.hxx"

    namespace JJCorrFitter
    {
        class CorrelationFunctionImpl
        {
            private:
                std::unique_ptr<InteractionTermImpl> m_psiSqrt;
                std::unique_ptr<SourceFunctionImpl> m_sourceFunc;

            public:
                CorrelationFunctionImpl(/* args */);
                virtual ~CorrelationFunctionImpl(){}
                virtual double Evaluate() = 0;
        };

    } // namespace JJCorrFitter
    

#endif