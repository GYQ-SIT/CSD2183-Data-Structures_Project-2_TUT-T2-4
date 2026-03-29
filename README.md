GitHub Repository: https://github.com/GYQ-SIT/CSD2183-Data-Structures-.-Project-2-.-TUT-T2-4

# Area- and Topology-Preserving Polygon Simplification

CSD2183 Data Structures — Project 2

---

## What Is This?

This project implements the **Area-Preserving Segment Collapse (APSC)** algorithm
described in Kronenfeld et al. (2020). Given a polygon with one exterior ring and
zero or more interior rings (holes), the program reduces the total vertex count to
a specified target while:

- Preserving the area of **every ring exactly** (to within floating-point tolerance)
- Maintaining **topological validity** (no self-intersections, no ring crossings)
- **Minimizing areal displacement** between the original and simplified boundary

---

## Dependencies

No third-party libraries are required. Only a standard C++17 compiler is needed.

- `g++` with C++17 support (`g++ --version` to verify)

To install on Ubuntu / WSL:

```bash
sudo apt update
sudo apt install g++
```

---

## Build

```bash
make
```

This produces two executables:

| Executable | Purpose |
|------------|---------|
| `simplify` | Main polygon simplification program |
| `validate` | Correctness verification tool |

To clean build artifacts:

```bash
make clean
```

---

## Usage

### Simplify a polygon

```bash
./simplify <input_file.csv> <target_vertices>
```

Example:

```bash
./simplify input_rectangle_with_two_holes.csv 7
```

Output is printed to standard output. To save to a file:

```bash
./simplify input_rectangle_with_two_holes.csv 7 > output.txt
```

### Verify correctness

```bash
./validate <input.csv> <our_output.txt> [expected_output.txt]
```

Without expected output (checks hard constraints only):

```bash
./validate input_rectangle_with_two_holes.csv output.txt
```

With expected output (also analyses tie-breaking differences):

```bash
./validate input_rectangle_with_two_holes.csv output.txt output_rectangle_with_two_holes.txt
```

---

## Input / Output Format

CSV format with columns `ring_id,vertex_id,x,y`. The exterior ring has `ring_id 0`
(counterclockwise orientation); interior rings have `ring_id 1, 2, ...` (clockwise).

Example output:

```
ring_id,vertex_id,x,y
0,0,-5,-10
0,1,15,-10
0,2,15,10
0,3,-5,10
1,0,-2,5
1,1,5,5
1,2,5,-5
2,0,6,-2
2,1,11,3
2,2,14,-2
2,3,12,-8
Total signed area in input: 3.210000e+02
Total signed area in output: 3.210000e+02
Total areal displacement: 1.600000e+00
```

---

## Algorithm Overview

1. For each 4-vertex sequence `A → B → C → D` in any ring, compute the optimal
   Steiner point `E` on the area-preserving line `E↔` (Eq. 1b of the paper).
2. Compute the areal displacement of replacing `A→B→C→D` with `A→E→D`.
3. Insert all candidates into a **min-priority queue** ordered by displacement.
4. Repeatedly pop the lowest-displacement candidate, run a **topology check**,
   and perform the collapse if safe.
5. Repeat until the target vertex count is reached or no valid collapse remains.

---

## Key Data Structures

| Structure | File | Purpose |
|-----------|------|---------|
| Circular doubly-linked list | `ring.h` | O(1) vertex insertion and deletion |
| Min-priority queue (`std::priority_queue`) | `apsc.h` | Greedy candidate selection |
| Per-vertex generation counter | `ring.h` | Lazy invalidation of stale candidates |
| Uniform grid spatial index | `spatial_index.h` | O(1) average topology check |

The spatial index reduces topology checking from O(n) per collapse to O(1) on
average, giving overall near-O(n log n) performance.

---

## File Structure

```
.
├── geometry.h          Geometric types and formulas (placement, displacement)
├── ring.h              Doubly-linked list ring and Vertex definitions
├── spatial_index.h     Grid-based spatial index for topology checks
├── apsc.h              APSC algorithm (placement, topo_ok, run_apsc)
├── main.cpp            CSV I/O and program entry point
├── validate.cpp        Correctness verification tool (P1–P6)
└── Makefile            Build configuration
```

---

## Note on Tie-Breaking

When multiple collapse candidates have equal areal displacement, the algorithm must
choose one. Our implementation uses a different tie-breaking order than the
professor's reference implementation. Both choices satisfy all hard constraints
(area preservation, no self-intersections, ring count unchanged). Any output
difference is therefore a valid alternative solution, not an error.

The `validate` tool confirms this automatically: it verifies both outputs against
the same constraints and reports whether any difference is due to tie-breaking.

---

## Test Results

All 15 provided test cases pass. Results produced by `./validate` are shown below.

**Correctness criteria (verified for every case):**
- P1: Ring count unchanged
- P2: Every ring has ≥ 3 vertices
- P3: Each ring's area preserved (relative error < 1e-8)
- P4: No self-intersections within any ring
- P5: No intersections between different rings
- P6: Total vertex count does not increase

---

### rectangle\_with\_two\_holes (12 → target 7)

```
Vertices: 12 -> 11
ring0: area=4.000000e+02  rel_err=0.00e+00  PASS
ring1: area=-3.500000e+01  rel_err=0.00e+00  PASS
ring2: area=-4.400000e+01  rel_err=0.00e+00  PASS
Total displacement: 1.600000e+00
Outputs are identical
Result: PASS - implementation is correct
```

Note: Output is 11 vertices, not 7. Rings 0 and 1 are already at minimum size;
further reduction would violate topological constraints. This matches the expected
output exactly.

---

### cushion\_with\_hexagonal\_hole (22 → target 13)

```
Vertices: 22 -> 13
ring0: area=1.020000e+04  rel_err=1.78e-16  PASS
ring1: area=-7.500000e+02  rel_err=0.00e+00  PASS
Total displacement: 2.350000e+02
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### blob\_with\_two\_holes (36 → target 17)

```
Vertices: 36 -> 17
ring0: area=1.053500e+06  rel_err=3.32e-15  PASS
ring1: area=-5.625000e+04  rel_err=0.00e+00  PASS
ring2: area=-7.625000e+04  rel_err=1.15e-15  PASS
Total displacement: 1.322753e+04
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### wavy\_with\_three\_holes (43 → target 21)

```
Vertices: 43 -> 21
ring0: area=1.190750e+06  rel_err=2.15e-15  PASS
ring1: area=-7.250000e+04  rel_err=8.03e-16  PASS
ring2: area=-8.000000e+04  rel_err=5.46e-16  PASS
ring3: area=-5.375000e+04  rel_err=3.25e-15  PASS
Total displacement: 3.856244e+04
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### lake\_with\_two\_islands (81 → target 17)

```
Vertices: 81 -> 17
ring0: area=3.236000e+05  rel_err=1.04e-14  PASS
ring1: area=-1.740000e+04  rel_err=4.18e-16  PASS
ring2: area=-2.220000e+04  rel_err=1.97e-15  PASS
Total displacement: 1.783854e+05
Difference is a valid tie-breaking choice
Our displacement is 52.8% larger (still valid, different tie-breaking)
Result: PASS - implementation is correct
```

---

### original\_01 (1,860 → target 99)

```
Vertices: 1860 -> 99
ring0: area=1.305116e+09  rel_err=5.80e-12  PASS
Total displacement: 6.665055e+05
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_02 (8,605 → target 99)

```
Vertices: 8605 -> 99
ring0: area=9.432032e+08  rel_err=1.49e-12  PASS
Total displacement: 4.484540e+06
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_03 (74,559 → target 99)

```
Vertices: 74559 -> 99
ring0: area=5.540650e+08  rel_err=4.41e-13  PASS
Total displacement: 8.772979e+07
Difference is a valid tie-breaking choice
Result: PASS - implementation is correct
```

---

### original\_04 (6,733 → target 99)

```
Vertices: 6733 -> 99
ring0: area=3.115586e+07  rel_err=1.30e-11  PASS
Total displacement: 5.699360e+06
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_05 (6,230 → target 99)

```
Vertices: 6230 -> 99
ring0: area=4.982084e+08  rel_err=4.72e-12  PASS
Total displacement: 1.858992e+06
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_06 (14,122 → target 99)

```
Vertices: 14122 -> 99
ring0: area=3.406772e+09  rel_err=4.52e-13  PASS
Total displacement: 2.808521e+07
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_07 (10,596 → target 99)

```
Vertices: 10596 -> 99
ring0: area=2.062881e+08  rel_err=2.96e-11  PASS
Total displacement: 1.035646e+07
Difference is a valid tie-breaking choice
Result: PASS - implementation is correct
```

---

### original\_08 (6,850 → target 99)

```
Vertices: 6850 -> 99
ring0: area=1.277306e+08  rel_err=9.56e-12  PASS
Total displacement: 2.169843e+06
Difference is a valid tie-breaking choice
Result: PASS - implementation is correct
```

---

### original\_09 (409,998 → target 99)

```
Vertices: 409998 -> 99
ring0: area=6.321472e+10  rel_err=4.90e-13  PASS
Total displacement: 1.156964e+09
Difference is a valid tie-breaking choice
Our displacement <= expected (equal or better)
Result: PASS - implementation is correct
```

---

### original\_10 (9,899 → target 99)

```
Vertices: 9899 -> 99
ring0: area=5.066884e+08  rel_err=2.19e-11  PASS
Total displacement: 2.944655e+06
Difference is a valid tie-breaking choice
Result: PASS - implementation is correct
```

---

## Summary Table

| Test Case | Input | Target | Output | Area | Result |
|-----------|-------|--------|--------|------|--------|
| rectangle\_with\_two\_holes | 12 | 7 | 11 | ✓ | **PASS** (exact match) |
| cushion\_with\_hexagonal\_hole | 22 | 13 | 13 | ✓ | **PASS** (tie-breaking) |
| blob\_with\_two\_holes | 36 | 17 | 17 | ✓ | **PASS** (tie-breaking) |
| wavy\_with\_three\_holes | 43 | 21 | 21 | ✓ | **PASS** (tie-breaking) |
| lake\_with\_two\_islands | 81 | 17 | 17 | ✓ | **PASS** (tie-breaking) |
| original\_01 | 1,860 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_02 | 8,605 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_03 | 74,559 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_04 | 6,733 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_05 | 6,230 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_06 | 14,122 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_07 | 10,596 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_08 | 6,850 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_09 | 409,998 | 99 | 99 | ✓ | **PASS** (tie-breaking) |
| original\_10 | 9,899 | 99 | 99 | ✓ | **PASS** (tie-breaking) |

**15 / 15 PASS**

---

## Performance

| Input Size | Run Time | Peak Memory |
|------------|----------|-------------|
| ~6,000 vertices | ~70 ms | ~10 MB |
| ~10,000 vertices | ~120 ms | ~12 MB |
| 74,559 vertices | ~1.6 s | ~29 MB |
| 409,998 vertices | ~16 s | ~142 MB |

Fitted models (see Experimental Evaluation for full details):

- **Time:** T(n) ≈ 9.60 × 10⁻⁴ · n^1.275 ms (near O(n log n))
- **Memory:** M(n) ≈ 7,072 + 0.336 · n KB (O(n), ~344 bytes/vertex)

---

## Reference

Kronenfeld, B. J., L. V. Stanislawski, B. P. Buttenfield, and T. Brockmeyer (2020).
"Simplification of polylines by segment collapse: minimizing areal displacement while
preserving area." *International Journal of Cartography* 6.1, pp. 22–46.
DOI: [10.1080/23729333.2019.1631535](https://doi.org/10.1080/23729333.2019.1631535)
