/**
 * @file InteractionTermSquareWell.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Squared modulus of the two particle wavefunction. Square well solution (quantum-statistics only).
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef InteractionTermSquareWell_hxx
    #define InteractionTermSquareWell_hxx

    #include "InteractionTermImpl.hxx"
    #include <cmath>
    #include "boost/math/special_functions/bessel.hpp"

    namespace JJCorrFitter
    {
        class InteractionTermSquareWell : public InteractionTermImpl
        {
            private:
                static constexpr double m_hBarC{197.3269602};
                double m_kStar;

            public:
                InteractionTermSquareWell(/* args */);
                ~InteractionTermSquareWell() = default;
                InteractionTermSquareWell(const InteractionTermSquareWell&) = delete;
                InteractionTermSquareWell& operator=(const InteractionTermSquareWell&) = delete;
                InteractionTermSquareWell(InteractionTermSquareWell&&) noexcept = default;
                InteractionTermSquareWell& operator=(InteractionTermSquareWell&&) noexcept = default;

                void SetParameters(const std::vector<double> &pars);
                void SetMomentum(double kStar) noexcept;
                void PopulateGrid() {}
                [[nodiscard]] double GetValue(double rStar, double cosTheta) noexcept;
        };

    } // namespace JJCorrFitter
    

#endif