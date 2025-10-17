#pragma once
#include <vector>
#include "Particle.hpp"

// Tiny interface to allow different backends (GLUT, SDL, Vulkan, etc.)
class Renderer
{
public:
    virtual ~Renderer() = default;
    virtual void draw(const std::vector<particle_t> &particles) = 0;
    virtual void processEvents() = 0;
    virtual bool shouldClose() const = 0;
};