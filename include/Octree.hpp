#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>
#include "Particle.hpp"

/// Axis-aligned cube (same half-size along all axes)
struct AABB
{
    vector3_t center;   // cube center
    float     half;     // half-size

    inline bool contains(const vector3_t& p) const
    {
        const float hx = half, hy = half, hz = half;
        return (p.x >= center.x - hx && p.x < center.x + hx) &&
               (p.y >= center.y - hy && p.y < center.y + hy) &&
               (p.z >= center.z - hz && p.z < center.z + hz);
    }
};

/// Lightweight Octree for Barnes–Hut (float-based).
/// - Leaf nodes store a *small bucket* of particle indices to avoid pathological
///   subdivision when many particles overlap numerically.
/// - Internal nodes aggregate total mass and center of mass (COM).
class Octree
{
public:
    // Tuneables
    struct BuildParams {
        int   bucket_size;
        int   max_depth;
        float bounds_pad;
        BuildParams() : bucket_size(8), max_depth(32), bounds_pad(1e-2f) {}
    };

    Octree() = default;

    // Build the tree from particles positions/masses
    void build(const std::vector<particle_t>& particles,
               const BuildParams& bp)
    {
        bp_    = bp;
        parts_ = &particles;

        // Compute a cubic root bounds around all particles
        AABB root = computeRootBounds_();
        nodes_.clear();
        nodes_.reserve(std::max<size_t>(particles.size(), 8u));
        root_ = newNode_(root);

        // Insert particles
        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            insert_(root_, i, 0);
        }
    }

    // Compute total gravitational acceleration on particle i using BH approx
    // theta: opening angle (0.3..1.0 typical). Smaller = more accurate.
    vector3_t accelerationOn(int i, float G, float eps2, float theta) const
    {
        vector3_t acc = {0, 0, 0};
        if (nodes_.empty()) return acc;
        traverseAccumulate_(root_, i, G, eps2, theta, acc);
        return acc;
    }

#ifdef USE_MPI
    // MPI-friendly flattened node for broadcast
    struct MpiTreeNode {
        float center[3];   // AABB center
        float half;        // AABB half size
        float com[3];      // center of mass
        float mass;        // total mass in node
        int   child[8];    // indices of children, -1 if none
        int   leafOffset;  // index into leafIndices array, -1 if not leaf
        int   leafCount;   // number of entries for this leaf
        int   isLeaf;      // 1 if leaf, 0 otherwise
    };

    // Flatten the tree into arrays suitable for MPI broadcast
    void exportToMpiTree(std::vector<MpiTreeNode>& outNodes,
                         std::vector<int>& outLeafIndices) const
    {
        outNodes.clear();
        outLeafIndices.clear();
        if (root_ < 0 || nodes_.empty()) return;

        outNodes.reserve(nodes_.size());
        std::vector<int> idMap(nodes_.size(), -1);
        exportNodeRecursive_(root_, outNodes, outLeafIndices, idMap);
    }
#endif

private:
    struct Node
    {
        AABB      box;
        float     mass = 0.0f;         // total mass in this node
        vector3_t com  = {0,0,0};      // center of mass
        bool      leaf = true;

        // Children by octant (bit-coded: x(1), y(2), z(4) relative to center)
        std::array<int,8> child;       // index into nodes_, -1 if absent

        // Leaf bucket: indices of particles in this cell
        std::vector<int> bucket;

        Node(const AABB& b) : box(b)
        {
            child.fill(-1);
        }
    };

    // Members
    std::vector<Node>              nodes_;
    int                            root_  = -1;
    const std::vector<particle_t>* parts_ = nullptr;
    BuildParams                    bp_;

    // Helpers
    int newNode_(const AABB& b)
    {
        nodes_.emplace_back(b);
        return static_cast<int>(nodes_.size()) - 1;
    }

    // Compute cubic bounds
    AABB computeRootBounds_() const
    {
        vector3_t mn = { +std::numeric_limits<float>::infinity(),
                         +std::numeric_limits<float>::infinity(),
                         +std::numeric_limits<float>::infinity() };
        vector3_t mx = { -std::numeric_limits<float>::infinity(),
                         -std::numeric_limits<float>::infinity(),
                         -std::numeric_limits<float>::infinity() };

        if (!parts_ || parts_->empty()) {
            return {{0,0,0}, 1.0f};
        }

        for (const auto& p : *parts_) {
            mn.x = std::min(mn.x, p.position.x);
            mn.y = std::min(mn.y, p.position.y);
            mn.z = std::min(mn.z, p.position.z);
            mx.x = std::max(mx.x, p.position.x);
            mx.y = std::max(mx.y, p.position.y);
            mx.z = std::max(mx.z, p.position.z);
        }

        vector3_t center = {(mn.x + mx.x) * 0.5f,
                            (mn.y + mx.y) * 0.5f,
                            (mn.z + mx.z) * 0.5f};

        float dx = mx.x - mn.x;
        float dy = mx.y - mn.y;
        float dz = mx.z - mn.z;
        float extent = std::max(dx, std::max(dy, dz));
        float half   = 0.5f * extent;

        // Handle degenerate cases (all particles same pos)
        if (!(half > 0.0f)) half = 1.0f;

        half *= (1.0f + bp_.bounds_pad);
        return {center, half};
    }

    inline int childIndex_(const AABB& b, const vector3_t& p) const
    {
        int idx = 0;
        if (p.x >= b.center.x) idx |= 1;
        if (p.y >= b.center.y) idx |= 2;
        if (p.z >= b.center.z) idx |= 4;
        return idx;
    }

    inline AABB childBox_(const AABB& b, int oct) const
    {
        const float h2 = b.half * 0.5f;
        vector3_t c = b.center;
        c.x += (oct & 1) ? +h2 : -h2;
        c.y += (oct & 2) ? +h2 : -h2;
        c.z += (oct & 4) ? +h2 : -h2;
        return {c, h2};
    }

    // Incremental COM update
    inline void accumulate_(Node& n, const particle_t& p)
    {
        const float newMass = n.mass + p.mass;
        if (newMass <= 0.0f) return;
        vector3_t num = n.com * n.mass + p.position * p.mass;
        n.com.x = num.x / newMass;
        n.com.y = num.y / newMass;
        n.com.z = num.z / newMass;
        n.mass  = newMass;
    }

    void insert_(int nodeId, int partIndex, int depth)
    {
        Node& n = nodes_[nodeId];
        const particle_t& p = (*parts_)[partIndex];

        // Expand mass / COM first
        accumulate_(n, p);

        if (n.leaf) {
            n.bucket.push_back(partIndex);
            // Subdivide if bucket too large and depth limit not reached
            if (static_cast<int>(n.bucket.size()) > bp_.bucket_size &&
                depth < bp_.max_depth)
            {
                subdivide_(nodeId, depth);
            }
            return;
        }

        // Internal node: forward to child
        const int ci = childIndex_(n.box, p.position);
        if (n.child[ci] == -1) {
            n.child[ci] = newNode_(childBox_(n.box, ci));
        }
        insert_(n.child[ci], partIndex, depth + 1);
    }

    void subdivide_(int nodeId, int depth)
    {
        Node& n = nodes_[nodeId];
        // Create children lazily as needed while redistributing
        const auto bucket = n.bucket; // copy: we will clear and re-insert
        n.bucket.clear();
        n.leaf = false;

        for (int pi : bucket)
        {
            const particle_t& p = (*parts_)[pi];
            const int ci = childIndex_(n.box, p.position);
            if (n.child[ci] == -1) {
                n.child[ci] = newNode_(childBox_(n.box, ci));
            }
            insert_(n.child[ci], pi, depth + 1);
        }
    }

    // Traversal for acceleration accumulation on particle i
    void traverseAccumulate_(int nodeId, int i, float G, float eps2,
                             float theta, vector3_t& acc) const
    {
        const Node& n = nodes_[nodeId];
        if (n.mass <= 0.0f) return;

        const vector3_t& pi = (*parts_)[i].position;

        // Vector from particle to node COM
        vector3_t r = { n.com.x - pi.x, n.com.y - pi.y, n.com.z - pi.z };
        float r2_true = r.x * r.x + r.y * r.y + r.z * r.z;

        // If leaf, sum direct interactions with contained particles (skip self)
        if (n.leaf) {
            for (int j : n.bucket) {
                if (j == i) continue;
                const particle_t& pj = (*parts_)[j];
                vector3_t rij = { pj.position.x - pi.x,
                                  pj.position.y - pi.y,
                                  pj.position.z - pi.z };
                float d2     = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z + eps2;
                float inv_r  = 1.0f / std::sqrt(d2);
                float inv_r3 = inv_r * inv_r * inv_r;
                float s      = G * pj.mass * inv_r3;
                acc.x += rij.x * s;
                acc.y += rij.y * s;
                acc.z += rij.z * s;
            }
            return;
        }

        // Opening criterion: if (s / d) < theta, use monopole approximation
        const float s = n.box.half * 2.0f;      // node width
        float d       = std::sqrt(r2_true);     // distance to COM
        if ((s / d) < theta) {
            float r2_soft = r2_true + eps2;
            float inv_r   = 1.0f / std::sqrt(r2_soft);
            float inv_r3  = inv_r * inv_r * inv_r;
            float sgm     = G * n.mass * inv_r3;
            acc.x += r.x * sgm;
            acc.y += r.y * sgm;
            acc.z += r.z * sgm;
            return;
        }

        // Otherwise recurse to children
        for (int c = 0; c < 8; ++c) {
            int ci = n.child[c];
            if (ci != -1) {
                traverseAccumulate_(ci, i, G, eps2, theta, acc);
            }
        }
    }

#ifdef USE_MPI
    void exportNodeRecursive_(int nodeId,
                              std::vector<MpiTreeNode>& outNodes,
                              std::vector<int>& outLeafIndices,
                              std::vector<int>& idMap) const
    {
        if (idMap[nodeId] != -1) return; // already exported

        const Node& n = nodes_[nodeId];

        MpiTreeNode m{};
        m.center[0] = n.box.center.x;
        m.center[1] = n.box.center.y;
        m.center[2] = n.box.center.z;
        m.half      = n.box.half;

        m.com[0]    = n.com.x;
        m.com[1]    = n.com.y;
        m.com[2]    = n.com.z;
        m.mass      = n.mass;

        m.leafOffset = -1;
        m.leafCount  = 0;
        m.isLeaf     = n.leaf ? 1 : 0;

        for (int i = 0; i < 8; ++i) {
            m.child[i] = -1;
        }

        int newIdx = static_cast<int>(outNodes.size());
        outNodes.push_back(m);
        idMap[nodeId] = newIdx;

        if (n.leaf)
        {
            if (!n.bucket.empty()) {
                m.leafOffset = static_cast<int>(outLeafIndices.size());
                m.leafCount  = static_cast<int>(n.bucket.size());
                for (int pi : n.bucket) {
                    outLeafIndices.push_back(pi);
                }
            }
            outNodes[newIdx] = m;
            return;
        }

        // Internal node: export children
        for (int oct = 0; oct < 8; ++oct) {
            int ci = n.child[oct];
            if (ci == -1) continue;
            exportNodeRecursive_(ci, outNodes, outLeafIndices, idMap);
            int childNewIdx = idMap[ci];
            outNodes[newIdx].child[oct] = childNewIdx;
        }
    }
#endif
};
