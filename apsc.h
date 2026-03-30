#pragma once
#include "geometry.h"
#include "ring.h"
#include "spatial_index.h"
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <optional>

// ============================================================
// Placement function (Paper Eq.1b + Gastner tie-breaking)
// Computes the optimal Steiner point E on the area-preserving
// line E<-> for the four-vertex sequence A, B, C, D.
// ============================================================
inline std::optional<Point> placement(const Point& A, const Point& B,
                                       const Point& C, const Point& D) {
    double a = D.y - A.y;
    double b = A.x - D.x;
    double c = -B.y*A.x + (A.y-C.y)*B.x + (B.y-D.y)*C.x + C.y*D.x;

    if (std::abs(a) < 1e-15 && std::abs(b) < 1e-15) return std::nullopt;

    double sB = signed_dist(B, A, D);
    double sC = signed_dist(C, A, D);
    double dB = std::abs(sB), dC = std::abs(sC);

    // Fig 6a: E-line coincides with AD -> place E at D (zero displacement)
    {
        double val_A = a*A.x + b*A.y + c;
        if ((dB < 1e-12 && dC < 1e-12) || val_A == 0.0)
            return D;
    }

    ELine L{a, b, c};
    bool use_AB;

    if (dB > 1e-12 && dC > 1e-12 && ((sB > 0) == (sC > 0))) {
        // B and C on the same side of AD (includes Fig 6c parallel case)
        // Gastner tie-breaking: use AB if B is not farther than C from AD
        use_AB = (dB <= dC);
    } else if (dB > 1e-12 && dC > 1e-12) {
        // B and C on opposite sides of AD
        // Use whichever side E-line is on (same side as B)
        Point ep = (std::abs(b) > std::abs(a))
                   ? Point{0.0, -c/b}
                   : Point{-c/a, 0.0};
        double sE = signed_dist(ep, A, D);
        use_AB = ((sB > 0) == (sE > 0));
    } else {
        // One endpoint on AD: compare distances
        use_AB = (dB < dC);
    }

    std::optional<Point> result;
    if (use_AB) {
        result = line_isect(L, A, B);
        if (!result) result = line_isect(L, C, D); // fallback if AB || E-line
    } else {
        result = line_isect(L, C, D);
        if (!result) result = line_isect(L, A, B); // fallback if CD || E-line
    }
    return result;
}

// ============================================================
// Areal displacement of collapse ABCD -> AED
// Handles the two cases where DE or EA crosses BC.
// ============================================================
inline double areal_displacement(const Point& A, const Point& B,
                                  const Point& C, const Point& D,
                                  const Point& E) {
    if (segs_cross_bool(D, E, B, C)) {
        Point Q = seg_isect_pt(D, E, B, C);
        const Point q1[4] = {A, B, Q, E};
        const Point q2[3] = {Q, C, D};
        return poly_abs_area(q1, 4) + poly_abs_area(q2, 3);
    }
    if (segs_cross_bool(E, A, B, C)) {
        Point P = seg_isect_pt(E, A, B, C);
        const Point q1[3] = {A, B, P};
        const Point q2[4] = {P, C, D, E};
        return poly_abs_area(q1, 3) + poly_abs_area(q2, 4);
    }
    return std::abs(tri_signed_area(E, B, C) + tri_signed_area(E, C, D));
}

// ============================================================
// Candidate: one potential collapse operation
// ============================================================
struct Candidate {
    double  disp;    // areal displacement of this collapse
    Vertex* B_vert;
    Vertex* C_vert;
    Ring*   ring;
    Point   E;       // position of the Steiner point
    int     gen;     // generation of B_vert when candidate was created

    // Min-heap ordering with deterministic tie-breaking
    bool operator>(const Candidate& o) const {
        if (std::abs(disp - o.disp) > 1e-10) return disp > o.disp;
        return B_vert > o.B_vert;
    }
};

// ============================================================
// Topology check using spatial index
//
// Checks that new edges AE and ED do not intersect any existing
// edge (other than the three edges A->B, B->C, C->D being deleted).
//
// Key fix: A_prev->A and D->D_next are NOT excluded.
// These edges survive the collapse and must be checked.
// segs_cross_bool only detects strict interior crossings,
// so shared endpoints A and D do not cause false positives,
// but a new edge that non-trivially crosses A_prev->A or D->D_next
// must be detected (this was the root cause of self-intersections
// in the blob test case).
// ============================================================
inline bool topo_ok(const Point& A, const Point& E, const Point& D,
                    SpatialIndex& idx,
                    Vertex* A_vert, Vertex* B_vert,
                    Vertex* C_vert, Vertex* D_vert) {
    auto id_pair = [](Vertex* u, Vertex* v) -> std::pair<Vertex*,Vertex*> {
        return u < v ? std::make_pair(u,v) : std::make_pair(v,u);
    };
    using VP = std::pair<Vertex*,Vertex*>;
    // Only exclude the three edges being deleted by this collapse
    const VP excl[3] = {
        id_pair(A_vert, B_vert),
        id_pair(B_vert, C_vert),
        id_pair(C_vert, D_vert),
    };
    auto is_excluded = [&](Vertex* u, Vertex* v) {
        VP p = id_pair(u, v);
        for (const VP& e : excl) if (e == p) return true;
        return false;
    };

    // Query only edges whose bbox overlaps the combined bbox of AE and ED
    BBox bbAE = seg_bbox(A, E);
    BBox bbED = seg_bbox(E, D);
    BBox bb = { std::min(bbAE.xmin, bbED.xmin), std::min(bbAE.ymin, bbED.ymin),
                std::max(bbAE.xmax, bbED.xmax), std::max(bbAE.ymax, bbED.ymax) };

    auto candidates = idx.query(bb);
    for (Vertex* v : candidates) {
        if (!v->alive || !v->next->alive) continue;
        if (is_excluded(v, v->next)) continue;
        const Point& p1 = v->pos;
        const Point& p2 = v->next->pos;
        if (segs_cross_bool(A, E, p1, p2)) return false;
        if (segs_cross_bool(E, D, p1, p2)) return false;
        // Check if an existing vertex lies strictly on a new edge
        if (point_on_seg(p1, A, E)) return false;
        if (point_on_seg(p1, E, D)) return false;
    }
    return true;
}

// ============================================================
// Main APSC algorithm
// ============================================================
struct APSCResult { double total_disp; };

inline APSCResult run_apsc(std::vector<Ring*>& all_rings, int target_n) {
    int total_verts = 0;
    for (Ring* r : all_rings) total_verts += r->size;
    if (total_verts <= target_n) return {0.0};

    // Build spatial index
    SpatialIndex idx;
    idx.rebuild(all_rings);

    using PQ = std::priority_queue<Candidate,
                                   std::vector<Candidate>,
                                   std::greater<Candidate>>;
    PQ heap;

    // Create a candidate for the sequence A -> B -> C -> D
    // and push it onto the heap. Increments B_vert->gen to
    // invalidate any previously existing candidate for the same B_vert.
    auto make_candidate = [&](Ring* ring, Vertex* A_vert) {
        if (ring->size < 4) return;
        Vertex* B_vert = A_vert->next;
        Vertex* C_vert = B_vert->next;
        Vertex* D_vert = C_vert->next;
        if (A_vert == C_vert || B_vert == D_vert || A_vert == D_vert) return;

        auto E_opt = placement(A_vert->pos, B_vert->pos,
                               C_vert->pos, D_vert->pos);
        if (!E_opt) return;

        double d = areal_displacement(A_vert->pos, B_vert->pos,
                                       C_vert->pos, D_vert->pos, *E_opt);
        ++B_vert->gen; // invalidate old candidates for this B_vert

        Candidate cand;
        cand.disp   = d;
        cand.B_vert = B_vert;
        cand.C_vert = C_vert;
        cand.ring   = ring;
        cand.E      = *E_opt;
        cand.gen    = B_vert->gen;
        heap.push(cand);
    };

    // Seed all initial candidates
    for (Ring* ring : all_rings)
        ring->for_each_vertex([&](Vertex* v) { make_candidate(ring, v); });

    double total_disp = 0.0;

    while (total_verts > target_n && !heap.empty()) {
        Candidate cand = heap.top(); heap.pop();

        // Lazy deletion: skip dead vertices or stale candidates
        if (!cand.B_vert->alive || !cand.C_vert->alive) continue;
        if (cand.gen != cand.B_vert->gen) continue;

        Vertex* B_vert = cand.B_vert;
        Vertex* C_vert = cand.C_vert;
        Vertex* A_vert = B_vert->prev;
        Vertex* D_vert = C_vert->next;
        Ring*   ring   = cand.ring;

        if (ring->size < 4) continue;

        const Point& A = A_vert->pos;
        const Point& D = D_vert->pos;
        const Point& E = cand.E;

        if (!topo_ok(A, E, D, idx, A_vert, B_vert, C_vert, D_vert))
            continue;

        // Update spatial index: remove edges that will change
        idx.remove_edge(A_vert->prev); // prev->A
        idx.remove_edge(A_vert);       // A->B
        idx.remove_edge(B_vert);       // B->C
        idx.remove_edge(C_vert);       // C->D

        // Perform the collapse
        total_disp += cand.disp;
        ring->do_collapse(B_vert, C_vert, E);
        Vertex* E_vert = B_vert; // B_vert now holds E's position
        --total_verts;

        // Re-insert updated edges into spatial index
        idx.insert_edge(A_vert->prev); // prev->A
        idx.insert_edge(A_vert);       // A->E
        idx.insert_edge(E_vert);       // E->D

        // Regenerate candidates for all affected 4-vertex sequences
        // Five start vertices cover all sequences that now include E_vert
        for (Vertex* start : {A_vert->prev->prev,
                               A_vert->prev,
                               A_vert,
                               E_vert,
                               D_vert}) {
            make_candidate(ring, start);
        }
    }

    return {total_disp};
}
