#pragma once

struct vector3_t
{
    float x;
    float y;
    float z;

    vector3_t operator-(const vector3_t &other) const
    {
        vector3_t result = {0};
        result.x = x - other.x;
        result.y = y - other.y;
        result.z = z - other.z;
        return result;
    }

    vector3_t &operator-=(const vector3_t &other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    vector3_t operator*(float multiplier) const
    {
        vector3_t result = {0};
        result.x = x * multiplier;
        result.y = y * multiplier;
        result.z = z * multiplier;
        return result;
    }

    vector3_t operator+(const vector3_t &other) const
    {
        vector3_t result = {0};
        result.x = x + other.x;
        result.y = y + other.y;
        result.z = z + other.z;
        return result;
    }

    vector3_t &operator+=(const vector3_t &other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
};

struct particle_t
{
    vector3_t position;
    vector3_t velocity;
    vector3_t acceleration;
    float mass;
};
/*
struct particleSoA
{
    std::vector<float> positions_x;
    std::vector<float> positions_y;
    std::vector<float> positions_z;
    vector3_t velocity;
    vector3_t acceleration;
    std::vector<float> masses;
};
*/