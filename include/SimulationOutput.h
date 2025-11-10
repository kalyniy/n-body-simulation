#pragma once
/*
struct SimulationStep
{
    float flat_positions[];
};
*/
struct SimulationOutputHeader
{
    size_t n_particles;
    size_t target_steps;
    size_t passed_steps;
};
/*
struct SimulationOutput
{
    SimulationOutputHeader header;
    SimulationStep steps[];
};
*/