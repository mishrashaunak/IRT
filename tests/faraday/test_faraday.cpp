#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "faraday.hpp"

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

    // Magnetic field B to be updated
    VecField<dim> B(layout, {Quantity::Bx, Quantity::By, Quantity::Bz});
    // Electric field E used in Faraday’s law
    VecField<dim> E(layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez});

    double dx = layout->cell_size(Direction::X);
    double dt = 0.01; // small timestep

    // Initialize E with known analytic functions
    // E.y = cos(x), E.z = sin(x)
    for (std::size_t i = layout->dual_dom_start(Direction::X);
         i <= layout->dual_dom_end(Direction::X); ++i)
    {
        double x = layout->coordinate(Direction::X, Quantity::Ey, i);
        E.y(i) = std::cos(x);
        E.z(i) = std::sin(x);

        // Start with zero B
        B.y(i) = 0.0;
        B.z(i) = 0.0;
    }

    Faraday<dim> faraday(layout, dt);
    // Update B using Faraday’s law: B^{n+1} = B^n - dt * curl(E)
    faraday(E, B);

    std::cout << "Testing Faraday's law implementation in 1D...\n\n";

    // Compare with analytical result: 
    // By^{n+1} = dt * dEz/dx, Bz^{n+1} = -dt * dEy/dx
    for (std::size_t i = layout->dual_dom_start(Direction::X);
         i < layout->dual_dom_end(Direction::X); ++i)
    {
        double x = layout->coordinate(Direction::X, Quantity::Ey, i);

        double dEz_dx = (std::sin(x + dx) - std::sin(x)) / dx;
        double dEy_dx = (std::cos(x + dx) - std::cos(x)) / dx;

        double expected_By = dt * dEz_dx;
        double expected_Bz = -dt * dEy_dx;

        std::cout << "i=" << i
                  << "\tB.y = " << B.y(i) << "\tExpected: " << expected_By
                  << "\tError: " << std::abs(B.y(i) - expected_By) << "\n";

        std::cout << "\ti=" << i
                  << "\tB.z = " << B.z(i) << "\tExpected: " << expected_Bz
                  << "\tError: " << std::abs(B.z(i) - expected_Bz) << "\n";
    }

    return 0;
}

