/**
 * @file LogLikelihoodTest.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Goodness-of-fit function for the minimizer which utilises the logarithmic maximum likelyhood method
 * @version 0.1.0
 * @date 2024-10-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LogLikelihoodTest_hxx
    #define LogLikelihoodTest_hxx

    #include "LikelihoodImpl.hxx"

    namespace JJCorrFitter
    {
        class LogLikelihoodTest : public LikelihoodImpl
        {
            private:
                std::unique_ptr<TH1> m_dataSignal, m_dataBackground;

                double CalculateLogLikelihood(const std::unique_ptr<TH1> &signal, const std::unique_ptr<TH1> &background, const std::unique_ptr<TH1> &model);

            public:
                LogLikelihoodTest(/* args */) = delete;
                LogLikelihoodTest(std::unique_ptr<TH1> &&data, std::unique_ptr<TH1> &&signal, std::unique_ptr<TH1> &&background, std::unique_ptr<CorrelationFunctionImpl> &&model);
                ~LogLikelihoodTest() = default;
                LogLikelihoodTest(const LogLikelihoodTest &) = delete;
                LogLikelihoodTest& operator=(const LogLikelihoodTest &) = delete;
                LogLikelihoodTest(LogLikelihoodTest &&) noexcept;
                LogLikelihoodTest& operator=(LogLikelihoodTest &&) noexcept;

                [[nodiscard]] std::function<double (const double *)> GetObjectiveFunction(const std::vector<std::size_t> &corrFuncIndexes,const std::vector<std::size_t> &srcIndexes,const std::vector<std::size_t> &psiIndexes);
                [[nodiscard]] std::size_t GetNParams() const;
        };

        inline std::size_t LogLikelihoodTest::GetNParams() const {return m_corrFunc->GetNParams();}
        
    } // namespace JJCorrFitter
    

#endif