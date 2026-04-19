#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>

/// Axis-aligned cube
struct AABB
{
    float cx, cy, cz;   // center
    float     half;     // half-size

    inline bool contains(float px, float py, float pz) const
    {
        return (px >= cx - half && px < cx + half) &&
               (py >= cy - half && py < cy + half) &&
               (pz >= cz - half && pz < cz + half);
    }
};

/// Lightweight Octree for Barnes–Hut.
/// Now takes SoA float pointers directly instead of particle_t.
class Octree
{
public:
    struct BuildParams {
        int   bucket_size = 8;
        int   max_depth   = 32;
        float bounds_pad  = 1e-2f;
    };

    // can be used by both GPU and MPI
    struct FlatTreeNode {
        float center[3];
        float half;
        float com[3];
        float mass;
        int   child[8];
        int   leafOffset;
        int   leafCount;
        int   isLeaf;
    };
    using MpiTreeNode = FlatTreeNode;

    struct FlatLeafParticle {
        float pos[3];
        float mass;
    };
    using MpiLeafParticle = FlatLeafParticle;

    Octree() = default;

    /// Build the tree from SoA position/mass arrays
    void build(const float* px, const float* py, const float* pz,
               const float* mass, int n, const BuildParams& bp)
    {
        bp_ = bp;
        px_ = px;
        py_ = py;
        pz_ = pz;
        mass_ = mass;
        n_ = n;

        AABB root = computeRootBounds_();
        nodes_.clear();
        nodes_.reserve(std::max(n, 8));
        root_ = newNode_(root);

        // Insert particles
        for (int i = 0; i < n; ++i) {
            insert_(root_, i, 0);
        }
    }

    // Compute total gravitational acceleration on particle i using BH approx
    // theta: opening angle (0.3..1.0 typical). Smaller = more accurate.
    void accelerationOn(int i, float G, float eps2, float theta,
                        float& ax, float& ay, float& az) const
    {
        ax = ay = az = 0.f;
        if (nodes_.empty()) return;
        traverseAccumulate_(root_, i, G, eps2, theta, ax, ay, az);
    }

    void exportToFlatTree(std::vector<FlatTreeNode>& outNodes,
                          std::vector<FlatLeafParticle>& outLeafParticles) const
    {
        outNodes.clear();
        outLeafParticles.clear();
        if (root_ < 0 || nodes_.empty()) return;

        outNodes.reserve(nodes_.size());
        std::vector<int> idMap(nodes_.size(), -1);
        exportNodeRecursive_(root_, outNodes, outLeafParticles, idMap);
    }

    void exportToMpiTree(std::vector<FlatTreeNode>& n,
                         std::vector<FlatLeafParticle>& l) const
    { exportToFlatTree(n, l); }

private:
    struct Node {
        AABB  box;
        float mass = 0.f;
        float com_x = 0.f, com_y = 0.f, com_z = 0.f;
        bool  leaf = true;
        std::array<int,8> child;
        std::vector<int> bucket;
        Node(const AABB& b) : box(b) { child.fill(-1); }
    };

    std::vector<Node> nodes_;
    int root_ = -1;
    const float* px_ = nullptr;
    const float* py_ = nullptr;
    const float* pz_ = nullptr;
    const float* mass_ = nullptr;
    int n_ = 0;
    BuildParams bp_;

    int newNode_(const AABB& b) {
        nodes_.emplace_back(b);
        return static_cast<int>(nodes_.size()) - 1;
    }

    AABB computeRootBounds_() const {
        if (n_ <= 0) return {0,0,0, 1.0f};
        float mnx = +1e30f, mny = +1e30f, mnz = +1e30f;
        float mxx = -1e30f, mxy = -1e30f, mxz = -1e30f;
        for (int i = 0; i < n_; ++i) {
            mnx = std::min(mnx, px_[i]); mxx = std::max(mxx, px_[i]);
            mny = std::min(mny, py_[i]); mxy = std::max(mxy, py_[i]);
            mnz = std::min(mnz, pz_[i]); mxz = std::max(mxz, pz_[i]);
        }
        float cx = (mnx+mxx)*0.5f, cy = (mny+mxy)*0.5f, cz = (mnz+mxz)*0.5f;
        float ext = std::max({mxx-mnx, mxy-mny, mxz-mnz});
        float h = ext * 0.5f;
        if (!(h > 0.f)) h = 1.0f;
        h *= (1.0f + bp_.bounds_pad);
        return {cx, cy, cz, h};
    }

    inline int childIndex_(const AABB& b, float x, float y, float z) const {
        int idx = 0;
        if (x >= b.cx) idx |= 1;
        if (y >= b.cy) idx |= 2;
        if (z >= b.cz) idx |= 4;
        return idx;
    }

    inline AABB childBox_(const AABB& b, int oct) const {
        float h2 = b.half * 0.5f;
        float cx = b.cx + ((oct & 1) ? +h2 : -h2);
        float cy = b.cy + ((oct & 2) ? +h2 : -h2);
        float cz = b.cz + ((oct & 4) ? +h2 : -h2);
        return {cx, cy, cz, h2};
    }

    inline void accumulate_(Node& n, int pi) {
        float pm = mass_[pi];
        float newMass = n.mass + pm;
        if (newMass <= 0.f) return;
        n.com_x = (n.com_x * n.mass + px_[pi] * pm) / newMass;
        n.com_y = (n.com_y * n.mass + py_[pi] * pm) / newMass;
        n.com_z = (n.com_z * n.mass + pz_[pi] * pm) / newMass;
        n.mass = newMass;
    }

    void insert_(int nodeId, int pi, int depth) {
        Node& n = nodes_[nodeId];
        accumulate_(n, pi);
        if (n.leaf) {
            n.bucket.push_back(pi);
            if (static_cast<int>(n.bucket.size()) > bp_.bucket_size && depth < bp_.max_depth)
                subdivide_(nodeId, depth);
            return;
        }
        int ci = childIndex_(n.box, px_[pi], py_[pi], pz_[pi]);
        if (n.child[ci] == -1) n.child[ci] = newNode_(childBox_(n.box, ci));
        insert_(n.child[ci], pi, depth + 1);
    }

    void subdivide_(int nodeId, int depth) {
        Node& n = nodes_[nodeId];
        auto bucket = n.bucket;
        n.bucket.clear();
        n.leaf = false;
        for (int pi : bucket) {
            int ci = childIndex_(n.box, px_[pi], py_[pi], pz_[pi]);
            if (n.child[ci] == -1) n.child[ci] = newNode_(childBox_(n.box, ci));
            insert_(n.child[ci], pi, depth + 1);
        }
    }

    void traverseAccumulate_(int nodeId, int i, float G, float eps2,
                             float theta, float& ax, float& ay, float& az) const
    {
        const Node& n = nodes_[nodeId];
        if (n.mass <= 0.f) return;
        float rx = n.com_x - px_[i], ry = n.com_y - py_[i], rz = n.com_z - pz_[i];
        float r2 = rx*rx + ry*ry + rz*rz;

        if (n.leaf) {
            for (int j : n.bucket) {
                if (j == i) continue;
                float dx = px_[j] - px_[i], dy = py_[j] - py_[i], dz = pz_[j] - pz_[i];
                float d2 = dx*dx + dy*dy + dz*dz + eps2;
                float inv = 1.f / std::sqrt(d2);
                float inv3 = inv * inv * inv;
                float s = G * mass_[j] * inv3;
                ax += dx * s; ay += dy * s; az += dz * s;
            }
            return;
        }

        float s_width = n.box.half * 2.f;
        float d = std::sqrt(r2);
        if ((s_width / d) < theta) {
            float r2s = r2 + eps2;
            float inv = 1.f / std::sqrt(r2s);
            float inv3 = inv * inv * inv;
            float s = G * n.mass * inv3;
            ax += rx * s; ay += ry * s; az += rz * s;
            return;
        }

        for (int c = 0; c < 8; ++c)
            if (n.child[c] != -1)
                traverseAccumulate_(n.child[c], i, G, eps2, theta, ax, ay, az);
    }

    void exportNodeRecursive_(int nodeId,
                              std::vector<FlatTreeNode>& outNodes,
                              std::vector<FlatLeafParticle>& outLeafParticles,
                              std::vector<int>& idMap) const
    {
        if (idMap[nodeId] != -1) return;
        const Node& n = nodes_[nodeId];
        FlatTreeNode m{};
        m.center[0] = n.box.cx;
        m.center[1] = n.box.cy;
        m.center[2] = n.box.cz;
        m.half = n.box.half;
        m.com[0] = n.com_x;
        m.com[1] = n.com_y;
        m.com[2] = n.com_z;
        m.mass = n.mass;
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
                m.leafOffset = static_cast<int>(outLeafParticles.size());
                m.leafCount  = static_cast<int>(n.bucket.size());

                for (int pi : n.bucket) {
                    FlatLeafParticle lp{};
                    lp.pos[0] = px_[pi];
                    lp.pos[1] = py_[pi];
                    lp.pos[2] = pz_[pi];
                    lp.mass   = mass_[pi];
                    outLeafParticles.push_back(lp);
                }
            }
            outNodes[newIdx] = m;
            return;
        }

        // Internal node: export children
        for (int oct = 0; oct < 8; ++oct) {
            int ci = n.child[oct];
            if (ci == -1) continue;
            exportNodeRecursive_(ci, outNodes, outLeafParticles, idMap);
            outNodes[newIdx].child[oct] = idMap[ci];
        }
    }
};
