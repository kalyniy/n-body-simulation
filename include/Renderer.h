#pragma once
#include "ParticleSoA.h"

class Renderer {
public:
    virtual ~Renderer() = default;
    virtual void draw(const ParticleSoA &particles) = 0;
    virtual void processEvents() = 0;
    virtual bool shouldClose() const = 0;
};