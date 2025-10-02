#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    Ampere(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
        auto const start = m_grid->dual_dom_start(Direction::X);
        auto const end   = m_grid->dual_dom_end(Direction::X) - 1;  // to avoid out-of-bounds at i+1

        for (int i = start; i <= end; ++i)
        {
            J.x(i) = 0.0;
            J.y(i) = -(B.z(i + 1) - B.z(i)) / dx;
            J.z(i) =  (B.y(i + 1) - B.y(i)) / dx;
        }

        // To handle boundary condn at last point
        J.x(end + 1) = 0.0;
        J.y(end + 1) = 0.0;
        J.z(end + 1) = 0.0;
    }
        
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP


