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
                int m_nqMax;
                double m_kStar;
                static constexpr int m_kStarMin{0};
                static constexpr int m_gevToMev{1000};
                static constexpr std::size_t m_nCosThetaPoints{30};
                static constexpr std::size_t m_nRPoints{61};
                
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