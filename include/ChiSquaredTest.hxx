/**
 * @file ChiSquaredTest.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Goodness-of-fit function for the minimizer whihc utilises the Chi-squared test
 * @version 0.2.0
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef ChiSquaredTest_hxx
    #define ChiSquaredTest_hxx

    #include "LikelihoodImpl.hxx"

    namespace JJCorrFitter
    {
        class ChiSquaredTest : public LikelihoodImpl
        {
            private:
                double CalculateChi2(const std::unique_ptr<TH1> &data, const std::unique_ptr<TH1> &model); // this is not ideal implementation; const unique_ptr does not guarantee const TH1! I should use shared_ptr or raw

            public:
                ChiSquaredTest(/* args */) = delete;
                ChiSquaredTest(std::unique_ptr<TH1> &&data, std::unique_ptr<CorrelationFunctionImpl> &&func);
                ~ChiSquaredTest() = default;
                ChiSquaredTest(const ChiSquaredTest&) = delete;
                ChiSquaredTest& operator=(const ChiSquaredTest&) = delete;
                ChiSquaredTest(ChiSquaredTest&&) noexcept;
                ChiSquaredTest& operator=(ChiSquaredTest&&) noexcept;

                [[nodiscard]] std::function<double (const double *)> GetObjectiveFunction(const std::vector<std::size_t> &corrFuncIndexes,const std::vector<std::size_t> &srcIndexes,const std::vector<std::size_t> &psiIndexes);
                [[nodiscard]] std::size_t GetNParams() const;
        };

        inline std::size_t ChiSquaredTest::GetNParams() const {return m_corrFunc->GetNParams();}

    } // namespace JJCorrFitter
    
#endif