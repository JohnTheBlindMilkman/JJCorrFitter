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
                static constexpr std::string_view m_coralParFile{"wfparameters.dat"};
                static constexpr int m_kStarMin{0};
                int m_nqMax, m_kStar;

            public:
                InteractionTermSchrodinger(/* args */);
                ~InteractionTermSchrodinger();
                InteractionTermSchrodinger(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger& operator=(const InteractionTermSchrodinger&) = delete;
                InteractionTermSchrodinger(InteractionTermSchrodinger&&) noexcept;
                InteractionTermSchrodinger& operator=(InteractionTermSchrodinger&&) noexcept;

                void SetParameters(int kStar) noexcept;
                [[nodiscard]] double GetValue(float rStar, float cosTheta);
        };

    } // namespace JJCorrFitter
    

#endif