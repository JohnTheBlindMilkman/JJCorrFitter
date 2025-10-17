/**
 * @file InteractionTermPhaseShift.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Squared modulus of the two particle wavefunction. Calculated using the phase-shift measurements done by Reid 1974???. Implemented in CorAL.
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef InteractionTermPhaseShift_hxx
    #define InteractionTermPhaseShift_hxx

    #include "InteractionTermImpl.hxx"

    #include "/home/jedkol/Downloads/CorAL/include/coral.h" // I will fix this later (said every programmer...)

    namespace JJCorrFitter
    {
        class InteractionTermPhaseShift : public InteractionTermImpl
        {
            private:
                CWaveFunction_pp_phaseshift m_waveFunction;
                const std::string m_coralParFile{"./wfparameters.dat"};
                static constexpr int m_kStarMin{0};
                static constexpr int m_gevToMev{1000};
                int m_nqMax;
                double m_kStar;

            public:
                InteractionTermPhaseShift(/* args */);
                ~InteractionTermPhaseShift();
                InteractionTermPhaseShift(const InteractionTermPhaseShift&) = delete;
                InteractionTermPhaseShift& operator=(const InteractionTermPhaseShift&) = delete;
                InteractionTermPhaseShift(InteractionTermPhaseShift&&) noexcept;
                InteractionTermPhaseShift& operator=(InteractionTermPhaseShift&&) noexcept;

                void SetParameters(const std::vector<double> &pars);
                void SetMomentum(double kStar);
                [[nodiscard]] double GetValue(double rStar, double cosTheta);
        };

    } // namespace JJCorrFitter
    

#endif