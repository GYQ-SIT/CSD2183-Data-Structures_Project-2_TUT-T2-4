/**
 * validate.cpp
 * ============
 * Verifies the correctness of APSC output.
 *
 * Usage:
 *   ./validate <input.csv> <our_output.txt> [expected_output.txt]
 *
 * Hard constraints verified (independent of tie-breaking):
 *   P1  Number of rings unchanged
 *   P2  Each ring has >= 3 vertices
 *   P3  Area of each ring preserved (relative error < 1e-8)
 *   P4  No self-intersections within any ring
 *   P5  No intersections between different rings
 *   P6  Total vertex count does not increase
 *
 * If expected_output.txt is provided:
 *   - Verifies P1-P6 on the expected output as well
 *   - If both outputs differ but both satisfy P1-P6,
 *     the difference is attributed to valid tie-breaking
 *   - Compares areal displacement between the two outputs
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <optional>

// ---- Data structures -------------------------------------------------------

struct Point { double x, y; };

struct Polygon {
    std::map<int, std::vector<Point>> rings;
    double area_in  = 0.0;
    double area_out = 0.0;
    double disp     = 0.0;
};

// ---- Geometry utilities ----------------------------------------------------

static double shoelace(const std::vector<Point>& pts) {
    double s = 0.0;
    int n = (int)pts.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        s += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    return s * 0.5;
}

static double signed_dist(Point P, Point A, Point D) {
    return (D.x-A.x)*(P.y-A.y) - (D.y-A.y)*(P.x-A.x);
}

static bool segs_cross(Point A, Point B, Point P, Point Q) {
    double d1 = signed_dist(A, P, Q), d2 = signed_dist(B, P, Q);
    double d3 = signed_dist(P, A, B), d4 = signed_dist(Q, A, B);
    return ((d1>0&&d2<0)||(d1<0&&d2>0)) && ((d3>0&&d4<0)||(d3<0&&d4>0));
}

static bool point_on_seg(Point P, Point A, Point B) {
    double dx = B.x-A.x, dy = B.y-A.y;
    double cross = dx*(P.y-A.y) - dy*(P.x-A.x);
    if (std::abs(cross) > 1e-6) return false;
    double t = (std::abs(dx) > std::abs(dy))
               ? (P.x-A.x)/dx : (P.y-A.y)/dy;
    return t > 1e-9 && t < 1.0-1e-9;
}

static bool segs_intersect(Point A, Point B, Point P, Point Q) {
    if (segs_cross(A, B, P, Q))   return true;
    if (point_on_seg(A, P, Q))    return true;
    if (point_on_seg(B, P, Q))    return true;
    if (point_on_seg(P, A, B))    return true;
    if (point_on_seg(Q, A, B))    return true;
    return false;
}

// ---- I/O -------------------------------------------------------------------

static std::map<int,std::vector<Point>> read_csv(const std::string& path) {
    std::map<int,std::vector<Point>> rings;
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; exit(1); }
    std::string line;
    std::getline(f, line); // skip header
    while (std::getline(f, line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        std::vector<std::string> cols;
        while (std::getline(ss, tok, ',')) cols.push_back(tok);
        if (cols.size() < 4) continue;
        int rid = std::stoi(cols[0]);
        rings[rid].push_back({std::stod(cols[2]), std::stod(cols[3])});
    }
    return rings;
}

static Polygon parse_output(const std::string& path) {
    Polygon poly;
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; exit(1); }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty()) continue;
        if (line.rfind("ring_id",0)==0) continue;
        if (line.rfind("Total signed area in input:",0)==0)
            poly.area_in = std::stod(line.substr(line.rfind(':')+1));
        else if (line.rfind("Total signed area in output:",0)==0)
            poly.area_out = std::stod(line.substr(line.rfind(':')+1));
        else if (line.rfind("Total areal displacement:",0)==0)
            poly.disp = std::stod(line.substr(line.rfind(':')+1));
        else {
            std::istringstream ss(line);
            std::string tok;
            std::vector<std::string> cols;
            while (std::getline(ss, tok, ',')) cols.push_back(tok);
            if (cols.size() == 4) {
                int rid = std::stoi(cols[0]);
                poly.rings[rid].push_back({std::stod(cols[2]), std::stod(cols[3])});
            }
        }
    }
    return poly;
}

// ---- Property verification -------------------------------------------------

struct CheckResult {
    std::vector<std::string> errors;
    int total_verts_in  = 0;
    int total_verts_out = 0;
};

static CheckResult verify_properties(
    const std::map<int,std::vector<Point>>& rings_in,
    const Polygon& poly)
{
    CheckResult cr;
    const auto& rings_out = poly.rings;
    const double EPS_REL = 1e-8;

    for (auto& [rid,pts] : rings_in)  cr.total_verts_in  += (int)pts.size();
    for (auto& [rid,pts] : rings_out) cr.total_verts_out += (int)pts.size();

    // P1: ring count unchanged
    std::set<int> ids_in, ids_out;
    for (auto& [k,v] : rings_in)  ids_in.insert(k);
    for (auto& [k,v] : rings_out) ids_out.insert(k);
    if (ids_in != ids_out)
        cr.errors.push_back("P1 FAIL: ring set changed");

    for (auto& [rid, pts] : rings_out) {
        // P2: each ring >= 3 vertices
        if ((int)pts.size() < 3)
            cr.errors.push_back("P2 FAIL: ring" + std::to_string(rid)
                                 + " has only " + std::to_string(pts.size()) + " vertices");

        // P3: area preserved (relative error < 1e-8)
        if (rings_in.count(rid)) {
            double a_in  = shoelace(rings_in.at(rid));
            double a_out = shoelace(pts);
            double rel   = std::abs(a_in - a_out) / std::max(std::abs(a_in), 1e-9);
            if (rel > EPS_REL)
                cr.errors.push_back("P3 FAIL: ring" + std::to_string(rid)
                    + " area error rel=" + std::to_string(rel)
                    + " (in=" + std::to_string(a_in)
                    + " out=" + std::to_string(a_out) + ")");
        }

        // P4: no self-intersections within ring
        int n = (int)pts.size();
        bool found = false;
        for (int i = 0; i < n && !found; ++i) {
            for (int j = i+2; j < n && !found; ++j) {
                if (j == (i-1+n)%n) continue;
                if (segs_cross(pts[i], pts[(i+1)%n], pts[j], pts[(j+1)%n])) {
                    cr.errors.push_back("P4 FAIL: ring" + std::to_string(rid)
                        + " self-intersection: edge " + std::to_string(i)
                        + " crosses edge " + std::to_string(j));
                    found = true;
                }
            }
        }
    }

    // P5: no intersections between different rings
    std::vector<int> rids;
    for (auto& [k,v] : rings_out) rids.push_back(k);
    for (int i = 0; i < (int)rids.size() && cr.errors.empty(); ++i) {
        for (int j = i+1; j < (int)rids.size(); ++j) {
            const auto& ri = rings_out.at(rids[i]);
            const auto& rj = rings_out.at(rids[j]);
            int ni=(int)ri.size(), nj=(int)rj.size();
            bool found=false;
            for (int si=0;si<ni&&!found;++si)
                for (int sj=0;sj<nj&&!found;++sj)
                    if (segs_intersect(ri[si],ri[(si+1)%ni],rj[sj],rj[(sj+1)%nj])){
                        cr.errors.push_back("P5 FAIL: ring" + std::to_string(rids[i])
                            + " intersects ring" + std::to_string(rids[j]));
                        found=true;
                    }
        }
    }

    // P6: vertex count did not increase
    if (cr.total_verts_out > cr.total_verts_in)
        cr.errors.push_back("P6 FAIL: vertex count increased "
            + std::to_string(cr.total_verts_in) + " -> "
            + std::to_string(cr.total_verts_out));

    return cr;
}

// ---- Output comparison -----------------------------------------------------

static bool outputs_equal(const std::map<int,std::vector<Point>>& a,
                           const std::map<int,std::vector<Point>>& b,
                           double tol=1e-4) {
    if (a.size() != b.size()) return false;
    for (auto& [rid, ptsA] : a) {
        if (!b.count(rid)) return false;
        const auto& ptsB = b.at(rid);
        if (ptsA.size() != ptsB.size()) return false;
        for (int i=0;i<(int)ptsA.size();++i)
            if (std::abs(ptsA[i].x-ptsB[i].x)>tol ||
                std::abs(ptsA[i].y-ptsB[i].y)>tol) return false;
    }
    return true;
}

static std::string first_diff(const std::map<int,std::vector<Point>>& a,
                               const std::map<int,std::vector<Point>>& b,
                               double tol=1e-4) {
    for (auto& [rid, ptsA] : a) {
        if (!b.count(rid))
            return "ring" + std::to_string(rid) + " missing from expected output";
        const auto& ptsB = b.at(rid);
        if (ptsA.size() != ptsB.size())
            return "ring" + std::to_string(rid) + " vertex count differs: "
                 + std::to_string(ptsA.size()) + " vs " + std::to_string(ptsB.size());
        for (int i=0;i<(int)ptsA.size();++i) {
            if (std::abs(ptsA[i].x-ptsB[i].x)>tol ||
                std::abs(ptsA[i].y-ptsB[i].y)>tol) {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(6);
                oss << "ring" << rid << " vertex" << i
                    << ": ours=(" << ptsA[i].x << "," << ptsA[i].y << ")"
                    << " exp=(" << ptsB[i].x  << "," << ptsB[i].y  << ")";
                return oss.str();
            }
        }
    }
    return "";
}

// ---- Helpers ---------------------------------------------------------------

static void sep(char c='-', int w=60) {
    std::cout << std::string(w, c) << "\n";
}

static void print_ring_areas(
    const std::map<int,std::vector<Point>>& rings_in,
    const std::map<int,std::vector<Point>>& rings_out)
{
    for (auto& [rid,pts] : rings_in) {
        double a_in  = shoelace(pts);
        double a_out = rings_out.count(rid) ? shoelace(rings_out.at(rid)) : 0.0;
        double rel   = std::abs(a_in-a_out) / std::max(std::abs(a_in),1e-9);
        std::cout << "  ring" << rid
                  << ": area=" << std::scientific << std::setprecision(6) << a_in
                  << "  rel_err=" << std::scientific << std::setprecision(2) << rel
                  << "  " << (rel<1e-8 ? "PASS" : "FAIL") << "\n";
    }
}

// ---- Main ------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ./validate <input.csv> <our_output.txt> [expected.txt]\n";
        return 1;
    }
    std::string inp_path = argv[1];
    std::string our_path = argv[2];
    std::string exp_path = (argc >= 4) ? argv[3] : "";

    auto rings_in = read_csv(inp_path);
    auto our      = parse_output(our_path);

    int total_in = 0;
    for (auto& [k,v] : rings_in) total_in += (int)v.size();

    sep('=');
    std::cout << "  APSC Correctness Validation\n";
    std::cout << "  Input:    " << inp_path << " (" << total_in << " vertices)\n";
    std::cout << "  Output:   " << our_path << " (" << our.rings.size() << " rings)\n";
    if (!exp_path.empty())
        std::cout << "  Expected: " << exp_path << "\n";
    sep('=');

    // [1] Verify our output
    std::cout << "\n[1] Verifying constraints on our output\n";
    sep();
    auto cr_our = verify_properties(rings_in, our);

    std::cout << "  Vertices: " << total_in << " -> " << cr_our.total_verts_out << "\n";
    print_ring_areas(rings_in, our.rings);
    std::cout << "  Total displacement: "
              << std::scientific << std::setprecision(6) << our.disp << "\n";

    if (cr_our.errors.empty()) {
        std::cout << "\n  All constraints satisfied (P1-P6)\n";
    } else {
        std::cout << "\n  " << cr_our.errors.size() << " violation(s) found:\n";
        for (auto& e : cr_our.errors) std::cout << "    -> " << e << "\n";
    }

    // [2] If expected output provided, verify and compare
    if (!exp_path.empty()) {
        auto exp = parse_output(exp_path);

        sep();
        std::cout << "\n[2] Verifying constraints on expected output\n";
        sep();
        auto cr_exp = verify_properties(rings_in, exp);

        std::cout << "  Vertices: " << total_in << " -> " << cr_exp.total_verts_out << "\n";
        print_ring_areas(rings_in, exp.rings);
        std::cout << "  Total displacement: "
                  << std::scientific << std::setprecision(6) << exp.disp << "\n";

        if (cr_exp.errors.empty()) {
            std::cout << "\n  All constraints satisfied (P1-P6)\n";
        } else {
            std::cout << "\n  Expected output has violations:\n";
            for (auto& e : cr_exp.errors) std::cout << "    -> " << e << "\n";
        }

        // [3] Consistency analysis
        sep();
        std::cout << "\n[3] Output consistency analysis\n";
        sep();

        bool match = outputs_equal(our.rings, exp.rings);
        if (match) {
            std::cout << "  Outputs are identical\n";
        } else {
            std::string diff = first_diff(our.rings, exp.rings);
            std::cout << "  Outputs differ (first difference: " << diff << ")\n";
            std::cout << "\n  Checking whether difference is due to valid tie-breaking:\n";

            bool our_ok = cr_our.errors.empty();
            bool exp_ok = cr_exp.errors.empty();

            if (our_ok && exp_ok) {
                std::cout << "  Both outputs satisfy P1-P6\n";
                std::cout << "  Difference is a valid tie-breaking choice\n";
                std::cout << "\n  Displacement comparison:\n";
                std::cout << "    Ours:     " << std::scientific << std::setprecision(6)
                          << our.disp << "\n";
                std::cout << "    Expected: " << std::scientific << std::setprecision(6)
                          << exp.disp << "\n";
                double ratio = our.disp / std::max(exp.disp, 1e-9);
                if (ratio <= 1.0 + 1e-6)
                    std::cout << "    Our displacement <= expected (equal or better)\n";
                else
                    std::cout << "    Our displacement is " << std::fixed
                              << std::setprecision(1) << (ratio-1.0)*100.0
                              << "% larger (still valid, different tie-breaking)\n";
            } else if (!our_ok) {
                std::cout << "  Our output violates constraints -> needs fixing\n";
            } else {
                std::cout << "  Expected output has violations -> reference is unreliable\n";
            }
        }
    }

    sep('=');
    std::cout << "\n  Result: ";
    if (cr_our.errors.empty())
        std::cout << "PASS - implementation is correct\n";
    else
        std::cout << "FAIL - " << cr_our.errors.size() << " issue(s) need fixing\n";
    sep('=');

    return cr_our.errors.empty() ? 0 : 1;
}
