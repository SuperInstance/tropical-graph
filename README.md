# tropical-graph

> **Graph operations in the tropical semiring. Max-plus algebra meets network theory.**

[![crates.io](https://img.shields.io/crates/v/tropical-graph.svg)](https://crates.io/crates/tropical-graph)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Tropical semiring (max-plus algebra) applied to graph operations. Tropical matrix powers compute shortest paths, tropical eigenvalues give max cycle mean, and tropical polyhedra define new graph-theoretic objects.

## The Tropical Semiring

In the **max-plus** (tropical) semiring:
- a ⊕ b = max(a, b)
- a ⊗ b = a + b

This swaps "addition" for "maximum" and "multiplication" for "addition". Why? Because:
- Matrix multiplication in tropical semiring = **shortest paths**
- Tropical matrix powers = **all-pairs shortest paths**
- Tropical eigenvalues = **maximum cycle mean**

## Quick Start

```rust
use tropical_graph::TropicalSemiring;

let ts = TropicalSemiring::new();
assert_eq!(ts.add(3.0, 5.0), 5.0);  // max
assert_eq!(ts.mul(3.0, 5.0), 8.0);  // sum
```

## Part of the SuperInstance math fleet

Tropical geometry appears across the fleet: `flux-algebra`, `tropical-neural`, `tropical-graph`. Together they form a complete tropical mathematics toolkit.

## License

MIT © [SuperInstance](https://github.com/SuperInstance)
