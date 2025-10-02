#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
    // TODO implement the Faraday class, hint - get inspiration from Ampere

public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension>& B)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            auto start = m_grid->dual_dom_start(Direction::X);
            auto end   = m_grid->dual_dom_end(Direction::X) - 1;

            for (int i = start; i <= end; ++i)
            {
                // Bx unchanged in 1D
                B.x(i) = 0.0;

                // By^{n+1} = By^n + dt * dEz/dx
                B.y(i) += m_dt * (E.z(i + 1) - E.z(i)) / dx;

                // Bz^{n+1} = Bz^n - dt * dEy/dx
                B.z(i) -= m_dt * (E.y(i + 1) - E.y(i)) / dx;
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;
};

#endif // HYBRIDIR_FARADAY_HPP

