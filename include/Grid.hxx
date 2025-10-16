/**
 * @file Grid.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Flat 3-dimensional grid. Storage is row-wise
 * @version 0.1.0
 * @date 2025-10-08
 * 
 * @copyright Copyright (c) 2025
 * 
 */

 #ifndef Grid_hxx
    #define Grid_hxx

    #include <vector>
    
    namespace JJCorrFitter
    {
        template <typename T>
        class Grid
        {
            private:
                std::size_t m_nX, m_nY, m_nZ;
                std::vector<T> m_data;

            public:
                Grid() : m_nX(0), m_nY(0), m_nZ(0), m_data({}) {}
                Grid(std::size_t nX, std::size_t nY, std::size_t nZ) : m_nX(nX), m_nY(nY), m_nZ(nZ), m_data(m_nX * m_nY * m_nZ, T()) {}
                [[nodiscard]] T& operator()(std::size_t x, std::size_t y, std::size_t z) noexcept {return m_data[x + m_nX * y + m_nX * m_nY * z];}
                [[nodiscard]] const T& at(std::size_t x, std::size_t y, std::size_t z) const 
                {
                    if (x >= m_nX || y >= m_nY || z >= m_nZ)
                    {
                        std::stringstream ss;
                        ss << "GRID<T>::operator(): index out of bounds: x = " << x << " (max is " << m_nX - 1 << "), y = " << y << " (max is " << m_nY - 1 << "), z = " << z << " (max is " << m_nZ - 1 << ")";
                        throw std::out_of_range(ss.str().c_str());
                    }

                    return m_data.at(x + m_nX * y + m_nX * m_nY * z);
                }
        };
    } // namespace JJCorrFitter
    
#endif