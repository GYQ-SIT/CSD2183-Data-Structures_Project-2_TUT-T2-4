#include "apsc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <memory>
#include <stdexcept>

// ============================================================
// CSV reader: parses ring_id,vertex_id,x,y format
// ============================================================
static std::map<int, std::vector<Point>> read_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);

    std::map<int, std::vector<Point>> rings;
    std::string line;
    std::getline(f, line); // skip header

    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> cols;
        while (std::getline(ss, tok, ',')) cols.push_back(tok);
        if (cols.size() < 4) continue;
        int rid = std::stoi(cols[0]);
        rings[rid].push_back({std::stod(cols[2]), std::stod(cols[3])});
    }
    return rings;
}

// ============================================================
// Coordinate formatter: print integers without decimal point,
// otherwise use enough precision to represent the value.
// ============================================================
static std::string fmt_coord(double v) {
    if (v == (long long)v && std::abs(v) < 1e14) {
        std::ostringstream oss;
        oss << (long long)v;
        return oss.str();
    }
    std::ostringstream oss;
    oss << std::setprecision(15) << v;
    std::string s = oss.str();
    if (s.find('.') != std::string::npos) {
        size_t last = s.find_last_not_of('0');
        if (last != std::string::npos && s[last] == '.')
            s = s.substr(0, last);
        else
            s = s.substr(0, last + 1);
    }
    return s;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    std::string input_path = argv[1];
    int target_n = std::stoi(argv[2]);

    std::map<int, std::vector<Point>> rings_pts;
    try {
        rings_pts = read_csv(input_path);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    if (rings_pts.empty()) {
        std::cerr << "No rings found in input\n";
        return 1;
    }

    // Compute total signed area of input
    double total_area_in = 0.0;
    for (auto& [rid, pts] : rings_pts) {
        double s = 0.0;
        int n = (int)pts.size();
        for (int i = 0; i < n; ++i)
            s += pts[i].x * pts[(i+1)%n].y - pts[(i+1)%n].x * pts[i].y;
        total_area_in += s * 0.5;
    }

    // Build Ring objects
    std::vector<std::unique_ptr<Ring>> ring_ptrs;
    std::vector<Ring*> all_rings;
    for (auto& [rid, pts] : rings_pts) {
        ring_ptrs.push_back(std::make_unique<Ring>(rid, pts));
        all_rings.push_back(ring_ptrs.back().get());
    }

    // Run APSC
    auto result = run_apsc(all_rings, target_n);

    // Compute total signed area of output
    double total_area_out = 0.0;
    for (Ring* r : all_rings)
        total_area_out += r->signed_area();

    // Print CSV output
    std::cout << "ring_id,vertex_id,x,y\n";
    for (Ring* r : all_rings) {
        int vid = 0;
        r->for_each_vertex([&](Vertex* v) {
            std::cout << r->ring_id << "," << vid++ << ","
                      << fmt_coord(v->pos.x) << ","
                      << fmt_coord(v->pos.y) << "\n";
        });
    }

    // Print statistics
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Total signed area in input: "  << total_area_in  << "\n";
    std::cout << "Total signed area in output: " << total_area_out << "\n";
    std::cout << "Total areal displacement: "    << result.total_disp << "\n";

    return 0;
}
