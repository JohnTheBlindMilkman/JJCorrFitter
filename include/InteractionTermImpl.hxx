/**
 * @file InteractionTermImpl.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Squared modulus of the two particle wavefunction implementation
 * @version 0.1.0
 * @date 2024-09-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef InteractionTermImpl_hxx
    #define InteractionTermImpl_hxx

    #include <memory>
    #include <string>
    #include <string_view>
    #include <vector>

    #include "Grid.hxx"

    namespace JJCorrFitter
    {
        class InteractionTermImpl
        {
            private:
                /* data */

            protected:
                std::size_t m_numberOfParams;
                std::string_view m_InteractionTermName;
                std::vector<double> m_qPoints, m_rPoints, m_ctPoints;
                Grid<double> m_grid;

            public:
                InteractionTermImpl(/* args */) = default;
                virtual ~InteractionTermImpl(){}
                InteractionTermImpl(const InteractionTermImpl&) = default;
                InteractionTermImpl& operator=(const InteractionTermImpl&) = default;
                InteractionTermImpl(InteractionTermImpl&&) noexcept = default;
                InteractionTermImpl& operator=(InteractionTermImpl&&) noexcept = default;

                virtual void SetParameters(const std::vector<double> &pars) = 0;
                virtual void SetMomentum(double kStar) = 0;
                void SetMomentumBins(std::vector<double> &&qBins) noexcept {m_qPoints = std::move(qBins);}
                void SetDistanceBins(std::vector<double> &&rBins) noexcept {m_rPoints = std::move(rBins);}
                void SetCosThetaBins(std::vector<double> &&ctBins) noexcept {m_ctPoints = std::move(ctBins);}
                virtual void PopulateGrid() = 0;
                [[nodiscard]] virtual double GetValue(double rStar, double cosTheta) = 0;
                [[nodiscard]] std::size_t GetNParams() const noexcept;
                [[nodiscard]] std::string_view GetInteractiontermName() const noexcept;
        };

        inline std::size_t InteractionTermImpl::GetNParams() const noexcept {return m_numberOfParams;}
        inline std::string_view InteractionTermImpl::GetInteractiontermName() const noexcept {return m_InteractionTermName;}

    } // namespace JJCorrFitter
    

#endif