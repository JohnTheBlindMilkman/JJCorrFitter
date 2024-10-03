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

    namespace JJCorrFitter
    {
        class InteractionTermSchrodinger : public InteractionTermImpl
        {
            private:
                std::unique_ptr<CWaveFunction> m_waveFunction;
                const std::string m_coralParFile{"./wfparameters.dat"};
                static constexpr int m_kStarMin{0};
                static constexpr int m_gevToMev{1000};
                int m_nqMax;
                float m_kStar;

            public:
                InteractionTermSchrodinger(/* args */);
                ~InteractionTermSchrodinger();
                InteractionTermSchrodinger(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger& operator=(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger(InteractionTermSchrodinger&&) noexcept;
                InteractionTermSchrodinger& operator=(InteractionTermSchrodinger&&) noexcept;

                void SetParameters(const std::vector<double> &pars);
                void SetMomentum(float kStar);
                [[nodiscard]] double GetValue(float rStar, float cosTheta);
        };

    } // namespace JJCorrFitter
    

#endif