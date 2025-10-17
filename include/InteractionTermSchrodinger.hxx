/**
 * @file InteractionTermSchrodinger.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Squared modulus of the two particle wavefunction. Calculated as a solution to the stationary Schoredinger equation with Reid potentnial. Implemented in CorAL.
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef InteractionTermSchrodinger_hxx
    #define InteractionTermSchrodinger_hxx

    #include "InteractionTermImpl.hxx"

    #include "/home/jedkol/Downloads/CorAL/include/coral.h" // I will fix this later (said every programmer...)

    #include "Grid.hxx"

    #include "boost/math/quadrature/gauss.hpp"
    #include "boost/math/quadrature/gauss_kronrod.hpp"

    namespace JJCorrFitter
    {
        class InteractionTermSchrodinger : public InteractionTermImpl
        {
            private:
                CWaveFunction_pp_schrod m_waveFunction;
                const std::string m_coralParFile{"./wfparameters.dat"};
                // Grid<double> m_grid;
                int m_nqMax;
                double m_kStar;
                // std::vector<float> m_nMomBins;
                static constexpr int m_kStarMin{0};
                static constexpr int m_gevToMev{1000};
                static constexpr std::size_t m_nCosThetaPoints{30};
                static constexpr std::size_t m_nRPoints{61};
                // static constexpr std::array<float,m_nCosThetaPoints> m_cosThetaPoints
                // {
                //     -0.996893484074649,
                //     -0.983668123279747,
                //     -0.960021864968307,
                //     -0.926200047429274,
                //     -0.882560535792052,
                //     -0.829565762382768,
                //     -0.767777432104826,
                //     -0.697850494793315,
                //     -0.620526182989242,
                //     -0.536624148142019,
                //     -0.447033769538089,
                //     -0.352704725530878,
                //     -0.254636926167889,
                //     -0.153869913608583,
                //     -0.0514718425553176,
                //     0.0514718425553176,
                //     0.153869913608583,
                //     0.254636926167889,
                //     0.352704725530878,
                //     0.447033769538089,
                //     0.536624148142019,
                //     0.620526182989242,
                //     0.697850494793315,
                //     0.767777432104826,
                //     0.829565762382768,
                //     0.882560535792052,
                //     0.926200047429274,
                //     0.960021864968307,
                //     0.983668123279747,
                //     0.996893484074649
                // };
                // static constexpr std::array<float,m_nRPoints> m_rPoints
                // {
                //     0.0128897487377344,
                //     0.077662898133763,
                //     0.209225078239885,
                //     0.408296918006321,
                //     0.672091937471844,
                //     0.99945337579231,
                //     1.39063888128599,
                //     1.84499881426814,
                //     2.3606673075023,
                //     2.93598660519868,
                //     3.56986916134847,
                //     4.26085594043079,
                //     5.00680410445402,
                //     5.80556419737934,
                //     6.65524843866933,
                //     7.5537376301671,
                //     8.49847339683432,
                //     9.48684542526892,
                //     10.5163691043409,
                //     11.5843962964495,
                //     12.6879883034555,
                //     13.8241557615477,
                //     14.9899686292401,
                //     16.182381861728,
                //     17.3981699431593,
                //     18.6340768458027,
                //     19.8868720829422,
                //     21.1532521597854,
                //     22.4298265508315,
                //     23.713203936117,
                //     25,
                //     26.2867960638829,
                //     27.5701734491684,
                //     28.8467478402145,
                //     30.1131279170577,
                //     31.3659231541972,
                //     32.6018300568406,
                //     33.8176181382719,
                //     35.0100313707598,
                //     36.1758442384522,
                //     37.3120116965444,
                //     38.4156037035504,
                //     39.483630895659,
                //     40.513154574731,
                //     41.5015266031656,
                //     42.4462623698328,
                //     43.3447515613306,
                //     44.1944358026206,
                //     44.9931958955459,
                //     45.7391440595692,
                //     46.4301308386515,
                //     47.0640133948013,
                //     47.6393326924976,
                //     48.1550011857318,
                //     48.609361118714,
                //     49.0005466242076,
                //     49.3279080625281,
                //     49.5917030819936,
                //     49.7907749217601,
                //     49.9223371018662,
                //     49.9871102512622
                // };

            public:
                InteractionTermSchrodinger(/* args */);
                ~InteractionTermSchrodinger();
                InteractionTermSchrodinger(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger& operator=(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger(InteractionTermSchrodinger&&) noexcept;
                InteractionTermSchrodinger& operator=(InteractionTermSchrodinger&&) noexcept;

                void SetParameters(const std::vector<double> &pars);
                void SetMomentum(double kStar);
                void PopulateGrid();
                [[nodiscard]] double GetValue(double rStar, double cosTheta);
        };

    } // namespace JJCorrFitter
    

#endif