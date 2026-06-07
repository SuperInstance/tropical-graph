# tropical-graph

> **Graph operations in the tropical semiring. Max-plus algebra meets network theory.**

[![crates.io](https://img.shields.io/crates/v/tropical-graph.svg)](https://crates.io/crates/tropical-graph)
[![docs.rs](https://docs.rs/tropical-graph/badge.svg)](https://docs.rs/tropical-graph)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust library implementing tropical semiring (max-plus algebra) operations for graph problems. Tropical matrix multiplication computes shortest paths, tropical matrix powers yield all-pairs shortest paths, tropical eigenvalues give the maximum cycle mean, and tropical polyhedra define new combinatorial objects. Turns the algebraic structure of the tropical semiring into efficient graph algorithms.

---

## Table of Contents

- [What is Tropical Algebra?](#what-is-tropical-algebra)
- [Why Does This Matter?](#why-does-this-matter)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
- [Mathematical Background](#mathematical-background)
- [Installation](#installation)
- [Related Crates](#related-crates)
- [License](#license)

---

## What is Tropical Algebra?

The **tropical semiring** (ℝ ∪ {−∞}, ⊕, ⊗) replaces ordinary addition and multiplication:

```
a ⊕ b = max(a, b)       ("tropical addition")
a ⊗ b = a + b           ("tropical multiplication")
```

This single swap — max for +, + for × — creates an entirely new algebraic universe where:

- **Matrix multiplication** = shortest path computation
- **Matrix powers** = all-pairs shortest paths (Floyd-Warshall)
- **Eigenvalues** = maximum cycle mean (Karp's algorithm)
- **Kleene star** = transitive closure (reachability)

```
Standard algebra:  (A × B)_{ij} = Σ_k A_{ik} × B_{kj}
Tropical algebra:  (A ⊗ B)_{ij} = min_k (A_{ik} + B_{kj})
                                       ↑ shortest path! ↑
```

The tropical semiring is not a curiosity — it's a fundamental algebraic structure that appears in optimization, scheduling, algebraic geometry, phylogenetics, and control theory.

## Why Does This Matter?

**For shortest paths**: Tropical matrix powers compute all-pairs shortest paths in O(n³) — the same as Floyd-Warshall, but expressed as simple matrix algebra. The algebraic formulation makes it composable and parallelizable.

**For scheduling**: The tropical semiring models project scheduling where the completion time of a task depends on the maximum (not sum) of its prerequisites. Critical path analysis is tropical algebra.

**For game theory**: Tropical games model turn-based games where each player minimizes the worst-case cost — naturally expressed in the min-plus semiring.

**For algebraic geometry**: Tropical geometry converts algebraic varieties into polyhedral complexes. The tropicalization of a curve is a metric graph — a deep connection between algebra and combinatorics.

**For phylogenetics**: The tropical semiring models evolutionary distance computation on phylogenetic trees, where branch lengths combine additively and distances take the minimum.

## Architecture

```
tropical-graph
│
├── TropicalSemiring           ← Core algebraic element
│   ├── new(v)                     Wrap a value
│   ├── add(&other)                max(self, other)
│   ├── mul(&other)                self + other
│   ├── pow(n)                     n-fold multiplication
│   ├── div(&other)                Subtraction (inverse of mul)
│   ├── zero()                     -∞ (additive identity)
│   └── one()                      0 (multiplicative identity)
│
├── TropicalMatrix             ← Matrices over the tropical semiring
│   ├── zeros(rows, cols)          All entries = -∞
│   ├── identity(n)                Diagonal = 0, off-diag = -∞
│   ├── matmul(&other)             Tropical matrix multiplication
│   ├── power(n)                   Tropical matrix power (shortest paths)
│   └── kleene_star()              A* = I ⊕ A ⊕ A² ⊕ ... (transitive closure)
│
├── TropicalShortestPath       ← Path computation
│   ├── apsp(adj)                  All-pairs shortest paths
│   ├── shortest_paths(adj)        Via Kleene star
│   └── sssp(adj, source)          Single-source shortest paths
│
├── TropicalEigenvalue         ← Spectral theory
│   └── max_cycle_mean(adj)        Karp's algorithm: max mean-weight cycle
│
└── TropicalPolyhedron         ← Tropical geometry
    ├── new()                      Empty constraint set
    ├── add_halfspace(coeffs, b)   Add constraint: max(coeffs) ≤ b
    ├── contains(point)            Feasibility check
    └── num_constraints()          Number of constraints
```

## Quick Start

```rust
use tropical_graph::{
    TropicalSemiring, TropicalMatrix,
    TropicalShortestPath, TropicalEigenvalue,
};

// Tropical arithmetic
let a = TropicalSemiring::new(3.0);
let b = TropicalSemiring::new(5.0);
assert_eq!(a.add(&b).0, 5.0);   // max(3, 5) = 5
assert_eq!(a.mul(&b).0, 8.0);   // 3 + 5 = 8

// Build an adjacency matrix (use -∞ for no edge)
let mut adj = TropicalMatrix::zeros(4, 4);
adj.set(0, 1, TropicalSemiring::new(2.0));  // 0→1 costs 2
adj.set(0, 2, TropicalSemiring::new(5.0));  // 0→2 costs 5
adj.set(1, 2, TropicalSemiring::new(1.0));  // 1→2 costs 1
adj.set(1, 3, TropicalSemiring::new(6.0));  // 1→3 costs 6
adj.set(2, 3, TropicalSemiring::new(3.0));  // 2→3 costs 3

// All-pairs shortest paths via tropical matrix power
let dist = TropicalShortestPath::apsp(&adj);
// dist[0][3] should be min(2+6, 5+3) = 8 (via 0→1→3 or 0→2→3)
// Actually 0→1→2→3 = 2+1+3 = 6 (shorter!)

// Single-source shortest paths from node 0
let from_zero = TropicalShortestPath::sssp(&adj, 0);
for (i, d) in from_zero.iter().enumerate() {
    println!("Distance 0→{}: {:?}", i, d.0);
}

// Maximum cycle mean (tropical eigenvalue)
let lambda = TropicalEigenvalue::max_cycle_mean(&adj);
println!("Max cycle mean: {:.4}", lambda);
```

## API Reference

### TropicalSemiring

| Method | Returns | Description |
|--------|---------|-------------|
| `new(v)` | `Self` | Wrap value in tropical semiring |
| `zero()` | `Self` | Additive identity (−∞) |
| `one()` | `Self` | Multiplicative identity (0) |
| `is_zero()` | `bool` | Is value −∞? |
| `add(&other)` | `Self` | max(self, other) |
| `mul(&other)` | `Self` | self + other |
| `pow(n)` | `Self` | n-fold tropical multiplication |
| `div(&other)` | `Option<Self>` | Subtraction (inverse of mul) |

### TropicalMatrix

| Method | Returns | Description |
|--------|---------|-------------|
| `zeros(rows, cols)` | `Self` | Matrix of −∞ |
| `identity(n)` | `Self` | Tropical identity (diag=0, rest=−∞) |
| `get(i, j)` | `TropicalSemiring` | Access entry |
| `set(i, j, v)` | `()` | Set entry |
| `matmul(&other)` | `Self` | Tropical matrix multiplication |
| `power(n)` | `Self` | n-th tropical power |
| `kleene_star()` | `Self` | Transitive closure |

### Shortest Paths

| Method | Returns | Description |
|--------|---------|-------------|
| `TropicalShortestPath::apsp(adj)` | `TropicalMatrix` | All-pairs shortest paths |
| `TropicalShortestPath::shortest_paths(adj)` | `TropicalMatrix` | Via Kleene star |
| `TropicalShortestPath::sssp(adj, src)` | `Vec<TropicalSemiring>` | Single-source shortest paths |

### Eigenvalues & Polyhedra

| Method | Returns | Description |
|--------|---------|-------------|
| `TropicalEigenvalue::max_cycle_mean(adj)` | `f64` | Maximum mean-weight cycle |
| `TropicalPolyhedron::new()` | `Self` | Empty constraint set |
| `poly.add_halfspace(coeffs, b)` | `()` | Add tropical inequality |
| `poly.contains(point)` | `bool` | Feasibility check |

## Mathematical Background

### The Tropical Semiring

The max-plus (tropical) semiring (ℝ ∪ {−∞}, ⊕, ⊗):

```
a ⊕ b = max(a, b)         Addition
a ⊗ b = a + b             Multiplication
ε = −∞                    Additive identity (zero element)
e = 0                     Multiplicative identity (unit element)
```

This forms a **commutative semiring** — all ring axioms hold except additive inverses (there's no subtraction in max).

### Tropical Matrix Algebra

Tropical matrix multiplication:

```
(A ⊗ B)_{ij} = min_k (A_{ik} + B_{kj})
```

The n-th tropical matrix power A^⊗n gives the shortest path of length exactly n between all pairs. The Kleene star A* = ⊕_{k=0}^{∞} A^⊗k gives shortest paths of any length.

### Tropical Eigenvalues

A tropical eigenvalue λ and eigenvector v satisfy:

```
A ⊗ v = λ ⊗ v
```

In conventional notation: for all i, min_j(A_{ij} + v_j) = λ + v_i.

**Karp's formula** gives the maximum cycle mean:

```
λ_max = max_i min_k (A^n_{ki} − A^{n-k}_{ii}) / k
```

This is the tropical analogue of the spectral radius.

### Tropical Polyhedra

A tropical halfspace is defined by:

```
max_i(a_i + x_i) ≤ b
```

The intersection of tropical halfspaces forms a tropical polyhedron — a piecewise-linear analogue of classical polyhedra. Tropical polyhedra appear in:
- Worst-case execution time analysis
- Mean payoff game strategies
- Phylogenetic tree spaces

## Installation

```bash
cargo add tropical-graph
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
tropical-graph = "0.1"
```

## Related Crates

Part of the **SuperInstance Exocortex** math fleet:

- **[sheaf-laplacian](https://github.com/SuperInstance/sheaf-laplacian)** — Sheaf Laplacian and Hodge decomposition
- **[graph-homology](https://github.com/SuperInstance/graph-homology)** — Clique complexes and Betti numbers
- **[cohomology-ring](https://github.com/SuperInstance/cohomology-ring)** — Cup products and cohomology operations
- **[persistent-agent](https://github.com/SuperInstance/persistent-agent)** — Persistent homology for agent behavior
- **[categorical-coordination](https://github.com/SuperInstance/categorical-coordination)** — Category theory for multi-agent systems

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project.
