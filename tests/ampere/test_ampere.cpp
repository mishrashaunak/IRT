#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "ampere.hpp"

#include <iostream>
#include <cmath>
#include <memory>
#include <cassert>

int main()
{
    constexpr std::size_t dim = 1;
    std::array<std::size_t, dim> grid_size = {100};
    std::array<double, dim> cell_size = {0.1};
    std::size_t nbr_ghosts = 1;

    auto layout = std::make_shared<GridLayout<dim>>(grid_size, cell_size, nbr_ghosts);

    VecField<dim> B(layout, {Quantity::Bx, Quantity::By, Quantity::Bz});
    VecField<dim> J(layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz});

    double dx = layout->cell_size(Direction::X);

    // Initialize B.y and B.z with known functions
    for (std::size_t i = layout->dual_dom_start(Direction::X);
         i <= layout->dual_dom_end(Direction::X); ++i)
    {
        double x = layout->coordinate(Direction::X, Quantity::By, i);
        B.y(i) = std::sin(x);
        B.z(i) = std::cos(x);
    }

    Ampere<dim> ampere(layout);
    ampere(B, J);

    std::cout << "Testing Ampere's law implementation in 1D...\n\n";

    for (std::size_t i = layout->dual_dom_start(Direction::X);
         i < layout->dual_dom_end(Direction::X); ++i)
    {
        double x = layout->coordinate(Direction::X, Quantity::By, i);

        // Central difference approximations of dBz/dx and dBy/dx
        double dBz_dx = (std::cos(x + dx) - std::cos(x)) / dx;
        double dBy_dx = (std::sin(x + dx) - std::sin(x)) / dx;

        std::cout << "i=" << i
                  << "\tJ.y = " << J.y(i) << "\tExpected: " << -dBz_dx
                  << "\tError: " << std::abs(J.y(i) + dBz_dx) << "\n";

        std::cout << "\ti=" << i
                  << "\tJ.z = " << J.z(i) << "\tExpected: " << dBy_dx
                  << "\tError: " << std::abs(J.z(i) - dBy_dx) << "\n";
    }

    return 0;
}

