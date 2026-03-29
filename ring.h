#pragma once
#include "geometry.h"
#include <vector>

// ============================================================
// Doubly-linked list representation of a polygon ring
// Supports O(1) collapse operations
// ============================================================

struct Vertex {
    int     id;       // original vertex id from input
    Point   pos;
    Vertex* prev;
    Vertex* next;
    bool    alive;
    int     gen;      // generation counter for stale candidate detection
    int     ring_id;

    Vertex(int id, double x, double y, int ring_id = 0)
        : id(id), pos(x, y), prev(nullptr), next(nullptr),
          alive(true), gen(0), ring_id(ring_id) {}
};

class Ring {
public:
    int     ring_id;
    int     size;     // number of active vertices
    Vertex* head;     // arbitrary active vertex (traversal start)

    std::vector<Vertex*> storage; // owns all vertex memory

    Ring(int ring_id, const std::vector<Point>& pts)
        : ring_id(ring_id), size((int)pts.size()), head(nullptr) {
        int n = (int)pts.size();
        storage.reserve(n);
        for (int i = 0; i < n; ++i)
            storage.push_back(new Vertex(i, pts[i].x, pts[i].y, ring_id));
        for (int i = 0; i < n; ++i) {
            storage[i]->prev = storage[(i - 1 + n) % n];
            storage[i]->next = storage[(i + 1) % n];
        }
        head = storage[0];
    }

    ~Ring() {
        for (auto v : storage) delete v;
    }

    Ring(const Ring&) = delete;
    Ring& operator=(const Ring&) = delete;

    // Iterate over all active vertices (circular)
    template<typename Func>
    void for_each_vertex(Func f) const {
        Vertex* v = head;
        do { f(v); v = v->next; } while (v != head);
    }

    // Collect all vertex positions
    std::vector<Point> pts_list() const {
        std::vector<Point> result;
        result.reserve(size);
        Vertex* v = head;
        do { result.push_back(v->pos); v = v->next; } while (v != head);
        return result;
    }

    // Signed area (shoelace)
    double signed_area() const {
        double s = 0.0;
        Vertex* v = head;
        do {
            s += v->pos.x * v->next->pos.y - v->next->pos.x * v->pos.y;
            v = v->next;
        } while (v != head);
        return s * 0.5;
    }

    // Collapse B_vert and C_vert into Steiner point E.
    // Reuses B_vert node (updates its position), removes C_vert from the list.
    // Returns the new E_vert (which is B_vert with updated position).
    Vertex* do_collapse(Vertex* B_vert, Vertex* C_vert, const Point& E) {
        Vertex* D_vert = C_vert->next;
        B_vert->pos = E;
        B_vert->next = D_vert;
        D_vert->prev = B_vert;
        C_vert->alive = false;
        --size;
        // Fix: update head if it pointed to the deleted vertex
        if (head == C_vert) head = D_vert;
        return B_vert;
    }
};
