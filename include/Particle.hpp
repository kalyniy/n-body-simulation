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

    vector3_t operator*(float multiplier) const
    {
        vector3_t result = {0};

        result.x = x * multiplier;
        result.y = y * multiplier;
        result.z = z * multiplier;

        return result;
    }
};

struct particle_t
{
    vector3_t position;
    float mass;
    float acceleration;
};