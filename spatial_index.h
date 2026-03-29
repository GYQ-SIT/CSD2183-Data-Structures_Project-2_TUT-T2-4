#pragma once
#include "geometry.h"
#include "ring.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>

// ============================================================
// Grid-based spatial index
// Stores edges (keyed by their 'from' Vertex*) in grid cells.
// Used to accelerate topo_ok: only query cells overlapping
// the bounding box of the new edges AE and ED.
// ============================================================

struct BBox {
    double xmin, ymin, xmax, ymax;
};

inline BBox seg_bbox(const Point& A, const Point& B) {
    return { std::min(A.x, B.x), std::min(A.y, B.y),
             std::max(A.x, B.x), std::max(A.y, B.y) };
}

inline bool bbox_overlap(const BBox& a, const BBox& b) {
    return a.xmax >= b.xmin && b.xmax >= a.xmin
        && a.ymax >= b.ymin && b.ymax >= a.ymin;
}

class SpatialIndex {
public:
    double cell_size = 1.0;
    double ox = 0.0, oy = 0.0;
    int    cols = 1, rows = 1;

    // cell key -> list of 'from' vertices (edge = v -> v->next)
    std::unordered_map<int, std::vector<Vertex*>> grid;

    void rebuild(const std::vector<Ring*>& rings) {
        grid.clear();
        double xmin = 1e18, ymin = 1e18, xmax = -1e18, ymax = -1e18;
        int n = 0;
        for (Ring* r : rings)
            r->for_each_vertex([&](Vertex* v) {
                xmin = std::min(xmin, v->pos.x);
                ymin = std::min(ymin, v->pos.y);
                xmax = std::max(xmax, v->pos.x);
                ymax = std::max(ymax, v->pos.y);
                ++n;
            });
        double dx = xmax - xmin, dy = ymax - ymin;
        cell_size = std::max((dx + dy) / std::max(1, (int)std::sqrt((double)n)), 1e-9);
        ox = xmin; oy = ymin;
        cols = std::max(1, (int)(dx / cell_size) + 2);
        rows = std::max(1, (int)(dy / cell_size) + 2);
        for (Ring* r : rings)
            r->for_each_vertex([&](Vertex* v) { insert_edge(v); });
    }

    void insert_edge(Vertex* v) {
        BBox bb = seg_bbox(v->pos, v->next->pos);
        for_cells(bb, [&](int key) { grid[key].push_back(v); });
    }

    void remove_edge(Vertex* v) {
        BBox bb = seg_bbox(v->pos, v->next->pos);
        for_cells(bb, [&](int key) {
            auto& vec = grid[key];
            vec.erase(std::remove(vec.begin(), vec.end(), v), vec.end());
        });
    }

    // Query all edges whose bounding box overlaps with bb
    std::vector<Vertex*> query(const BBox& bb) const {
        std::vector<Vertex*> result;
        int c0 = std::max(0, (int)((bb.xmin - ox) / cell_size));
        int c1 = std::min(cols - 1, (int)((bb.xmax - ox) / cell_size) + 1);
        int r0 = std::max(0, (int)((bb.ymin - oy) / cell_size));
        int r1 = std::min(rows - 1, (int)((bb.ymax - oy) / cell_size) + 1);
        for (int r = r0; r <= r1; ++r)
            for (int c = c0; c <= c1; ++c) {
                auto it = grid.find(r * cols + c);
                if (it != grid.end())
                    for (Vertex* v : it->second)
                        result.push_back(v);
            }
        return result;
    }

private:
    template<typename F>
    void for_cells(const BBox& bb, F f) {
        int c0 = std::max(0, (int)((bb.xmin - ox) / cell_size));
        int c1 = std::min(cols - 1, (int)((bb.xmax - ox) / cell_size) + 1);
        int r0 = std::max(0, (int)((bb.ymin - oy) / cell_size));
        int r1 = std::min(rows - 1, (int)((bb.ymax - oy) / cell_size) + 1);
        for (int r = r0; r <= r1; ++r)
            for (int c = c0; c <= c1; ++c)
                f(r * cols + c);
    }
};
