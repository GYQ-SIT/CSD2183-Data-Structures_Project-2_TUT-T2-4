#pragma once
#include <cmath>
#include <optional>

// ============================================================
// Basic geometric types and functions
// All formulas verified against Python prototype
// ============================================================

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
};

// Signed area of triangle (positive = CCW)
inline double triangle_area(const Point& i, const Point& j, const Point& k) {
    return 0.5 * ((j.x - i.x) * (k.y - i.y) - (k.x - i.x) * (j.y - i.y));
}

// Signed distance of P from directed line A->D (positive = left side)
inline double signed_dist(const Point& P, const Point& A, const Point& D) {
    return (D.x - A.x) * (P.y - A.y) - (D.y - A.y) * (P.x - A.x);
}

// Coefficients of the area-preserving line E<-> (Eq. 1b in Kronenfeld et al.)
struct ELine { double a, b, c; };

inline ELine compute_E_line(const Point& A, const Point& B,
                             const Point& C, const Point& D) {
    double a = D.y - A.y;
    double b = A.x - D.x;
    double c = -B.y * A.x + (A.y - C.y) * B.x
               + (B.y - D.y) * C.x + C.y * D.x;
    return {a, b, c};
}

// Intersection of line ax+by+c=0 with the line through P1 and P2 (no range check)
inline std::optional<Point> line_isect(const ELine& L,
                                        const Point& P1, const Point& P2) {
    double dx = P2.x - P1.x, dy = P2.y - P1.y;
    double d = L.a * dx + L.b * dy;
    if (std::abs(d) < 1e-12) return std::nullopt;
    double t = -(L.a * P1.x + L.b * P1.y + L.c) / d;
    return Point{P1.x + t * dx, P1.y + t * dy};
}

// Strict interior intersection of segments P1P2 and P3P4 (excluding shared endpoints)
struct CrossResult { double t, s; Point pt; };

inline std::optional<CrossResult> seg_cross(const Point& P1, const Point& P2,
                                             const Point& P3, const Point& P4) {
    double dx1 = P2.x - P1.x, dy1 = P2.y - P1.y;
    double dx2 = P4.x - P3.x, dy2 = P4.y - P3.y;
    double denom = dx2 * dy1 - dx1 * dy2;
    if (std::abs(denom) < 1e-12) return std::nullopt;
    double t = (dx2 * (P3.y - P1.y) - dy2 * (P3.x - P1.x)) / denom;
    double s = (dx1 * (P3.y - P1.y) - dy1 * (P3.x - P1.x)) / denom;
    constexpr double EPS = 1e-9;
    if (t > EPS && t < 1.0 - EPS && s > EPS && s < 1.0 - EPS) {
        return CrossResult{t, s, {P1.x + t * dx1, P1.y + t * dy1}};
    }
    return std::nullopt;
}

// Returns true if segments AB and PQ cross strictly (no shared endpoints)
inline bool segs_cross_bool(const Point& A, const Point& B,
                              const Point& P, const Point& Q) {
    double d1 = signed_dist(A, P, Q), d2 = signed_dist(B, P, Q);
    double d3 = signed_dist(P, A, B), d4 = signed_dist(Q, A, B);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0))
        && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

// Intersection point of AB and PQ (assumes segs_cross_bool is true)
inline Point seg_isect_pt(const Point& A, const Point& B,
                            const Point& P, const Point& Q) {
    double d1 = signed_dist(A, P, Q), d2 = signed_dist(B, P, Q);
    double t = d1 / (d1 - d2);
    return {A.x + t * (B.x - A.x), A.y + t * (B.y - A.y)};
}

// Returns true if point P lies strictly inside segment AB (not at endpoints)
inline bool point_on_seg(const Point& P, const Point& A, const Point& B) {
    double dx = B.x - A.x, dy = B.y - A.y;
    double cross = dx * (P.y - A.y) - dy * (P.x - A.x);
    if (std::abs(cross) > 1e-9) return false;
    double t;
    if (std::abs(dx) > std::abs(dy))
        t = (std::abs(dx) > 1e-12) ? (P.x - A.x) / dx : 0.0;
    else
        t = (std::abs(dy) > 1e-12) ? (P.y - A.y) / dy : 0.0;
    return t > 1e-9 && t < 1.0 - 1e-9;
}

// Absolute area of polygon with n vertices (shoelace formula)
inline double poly_abs_area(const Point* pts, int n) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        s += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    return std::abs(s) * 0.5;
}

// Signed area of triangle EBC (positive = CCW)
inline double tri_signed_area(const Point& E, const Point& B, const Point& C) {
    return 0.5 * ((B.x - E.x) * (C.y - E.y) - (C.x - E.x) * (B.y - E.y));
}
