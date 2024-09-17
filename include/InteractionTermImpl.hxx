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

    namespace JJCorrFitter
    {
        class InteractionTermImpl
        {
            private:
                /* data */

            protected:
                float m_kStar,m_cosTheta;

            public:
                InteractionTermImpl(/* args */) : m_kStar(1.f), m_cosTheta(0.f) {}
                virtual ~InteractionTermImpl(){}
                InteractionTermImpl(const InteractionTermImpl&) = default;
                InteractionTermImpl& operator=(const InteractionTermImpl&) = default;
                InteractionTermImpl(InteractionTermImpl&&) noexcept = default;
                InteractionTermImpl& operator=(InteractionTermImpl&&) noexcept = default;

                virtual void SetParameters(float kStar) noexcept = 0;
                [[nodiscard]] virtual double GetValue(float rStar, float cosTheta) = 0;
        };

    } // namespace JJCorrFitter
    

#endif