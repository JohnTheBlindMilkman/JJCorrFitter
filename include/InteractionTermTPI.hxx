/**
 * @file InteractionTermTPI.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Squared modulus of the two particle wavefunction. Calculated with Lednicky-Lyuboshitz approach. Implemented in HAL and TPI.
 * @version 0.1.0
 * @date 2024-09-27
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef InteractionTermTPI_hxx
    #define InteractionTermTPI_hxx

    #include <complex>
    #include <iostream>
    #include <array>

    #include "Math/Math.h"

    #include "InteractionTermImpl.hxx"

    namespace JJCorrFitter
    {
        class InteractionTermTPI : public InteractionTermImpl
        {
            public:
                enum class SpinState{None,SpinAveraged,SpinSeparated};

            private:
                /**
                 * Calculates the confluent hypergeometric function F from single orientation of cos(theta*) for non-symmetrized wave-function
                 * (non-identical particles)
                 * @param ffp
                 * @param sign
                 */
                constexpr void GetFFsingle(float rStar, float cosTheta, std::complex<long double> &ffp, int sign = 1) const;
                /**
                 * Calculates the confluent hypergeometric function for two orientations of cos(theta*) for symmetrized wave-function
                 * (identical particles)
                 * @param ffp
                 * @param ffm
                 */
                constexpr void GetFFdouble(float rStar, float cosTheta, std::complex<long double> &ffp, std::complex<long double> &ffm) const;
                /**
                 * Calculates G~ function
                 * @param eta
                 * @param rho
                 * @param hfun
                 * @return
                 */
                constexpr std::complex<long double> GetG(long double eta, long double rho, long double hfun) const;
                constexpr long double Chiim(long double eta) const;
                /**
                 * Calculates H function for strong interaction
                 * @param eta
                 * @return
                 */
                constexpr long double GetH(long double eta) const;
                constexpr void Getfc(long double kstar, long double eta, long double hfun, std::complex<long double> &fcs, std::complex<long double> &fct) const;  // TODO
                constexpr std::pair<long double,long double> Bfunpfun(long double eta, long double rho) const;
                constexpr double Funeh(double xarg, double rad, double alfa) const;
                constexpr double Funex(double xarg, double rad) const;
                constexpr void InitializeGamow();
                constexpr double Gamow(double arg) const;
                constexpr double GetQuantumCoulombStrong(float rStar, float cosTheta);
                /**
                 * @brief Calculates weight for identical bosons
                 * 
                 * @return double
                 */
                constexpr double GetQuantumCoulomb(float rStar, float cosTheta);

                static constexpr float m_gevToFm{5.0677302};
                static constexpr float m_fmToGeV{1./m_gevToFm};
                static constexpr float m_mevToGev{0.001};
                static constexpr double m_pi{3.141592653589793238};
                static constexpr int fCoulombSteps{170};

                SpinState m_spinState;
                std::complex<long double> fD0s, fF0s, fD0t, fF0t;
                long double fPionac, fOneoveracsq, fTwopioverac, fCoulqscpart, fEuler, fF0, fD0;
                int fTwospin, fWritegrps, fPcount;
                double m_kStar, fRStarOutS, fRStarSideS, fRStarLongS, fRStarS, fRStarOut, fRStarSide, fRStarLong, fRStar, fKStarOut, fKStarSide, fKStarLong, fKStar;

            public:
                InteractionTermTPI() = delete;
                explicit InteractionTermTPI(SpinState state);
                ~InteractionTermTPI() = default;
                InteractionTermTPI(const InteractionTermTPI &) = delete;
                InteractionTermTPI& operator=(const InteractionTermTPI &) = delete;
                InteractionTermTPI(InteractionTermTPI &&) noexcept = default;
                InteractionTermTPI& operator=(InteractionTermTPI &&) noexcept = default;

                void SetParameters(const std::vector<double> &pars);
                void SetMomentum(float kStar);
                [[nodiscard]] double GetValue(float rStar, float cosTheta);
        };
        
        inline constexpr double InteractionTermTPI::Gamow(double arg) const 
        {
            long double eta = fTwopioverac / arg;
            return (eta) *1.0 / (exp(eta) - 1.0);
        }
        inline constexpr long double InteractionTermTPI::Chiim(long double eta) const { return Gamow(1.0 / (eta * fPionac)) / (2.0 * eta); }
        inline constexpr double InteractionTermTPI::Funeh(double xarg, double rad, double alfa) const { return exp(-sqrt(xarg * xarg / (rad * rad) + alfa * alfa)); }
        inline constexpr double InteractionTermTPI::Funex(double xarg, double rad) const { return exp(-xarg / rad); }

    } // namespace JJCorrFitter
    
#endif