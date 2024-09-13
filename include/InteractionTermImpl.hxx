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
            public:
                InteractionTermImpl(/* args */);
                virtual ~InteractionTermImpl(){}
                InteractionTermImpl(const InteractionTermImpl&) = delete;
                InteractionTermImpl& operator=(const InteractionTermImpl&) = delete;
                InteractionTermImpl(InteractionTermImpl&&) noexcept = default;
                InteractionTermImpl& operator=(InteractionTermImpl&&) noexcept = default;
                [[nodiscard]] virtual double GetValue(int kStar, float rStar, float cosTheta) = 0;
        };

    } // namespace JJCorrFitter
    

#endif