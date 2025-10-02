#include "population.hpp"
#include "moments.hpp"

#include <iostream>

int main()
{
    constexpr std::size_t dim = 1;
    std::array<std::size_t, dim> grid_size = {10};
    std::array<double, dim> cell_size = {1.0};
    std::size_t nbr_ghosts = 1;

    auto layout = std::make_shared<GridLayout<dim>>(grid_size, cell_size, nbr_ghosts);

    // Create a test population
    Population<dim> pop("test", layout);

    // Load 1 particle at cell center with fixed velocity
    Particle<dim> p;
    p.position[0] = layout->cell_coordinate(Direction::X, 5); // middle of grid
    p.v = {1.0, 0.0, 0.0};
    p.mass = 1.0;
    p.charge = 1.0;
    p.weight = 2.0;

    pop.particles().push_back(p);

    // Deposit moments
    pop.deposit();

    // Print density and flux
    std::cout << "Density:\n";
    for (std::size_t i = 0; i < pop.density().data().size(); ++i)
        std::cout << "n[" << i << "] = " << pop.density()(i) << "\n";

    std::cout << "\nFlux (Vx):\n";
    for (std::size_t i = 0; i < pop.flux().x.data().size(); ++i)
        std::cout << "Vx[" << i << "] = " << pop.flux().x(i) << "\n";


    // Wrap in vector for total_density and bulk_velocity
    std::vector<Population<dim>> pops = {pop};

    Field<dim> N(layout->allocate(Quantity::N), Quantity::N);
    total_density(pops, N);

    VecField<dim> V(layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz});
    bulk_velocity(pops, N, V);

    std::cout << "\nTotal Density N:\n";
    for (std::size_t i = 0; i < N.data().size(); ++i)
        std::cout << "N[" << i << "] = " << N(i) << "\n";

    std::cout << "\nBulk Velocity Vx:\n";
    for (std::size_t i = 0; i < V.x.data().size(); ++i)
        std::cout << "Vx[" << i << "] = " << V.x(i) << "\n";

    return 0;
}

