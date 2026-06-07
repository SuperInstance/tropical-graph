//! # tropical-graph
//!
//! Tropical semiring (max-plus) algebra for graph operations including
//! shortest paths via tropical matrix powers, tropical eigenvalues
//! (maximum cycle mean), and tropical polyhedra.

use std::fmt;

/// Tropical infinity (represents -∞ for max-plus).
pub const TROPICAL_NEG_INF: f64 = f64::NEG_INFINITY;
/// Tropical zero element for multiplication (additive identity in max-plus).
pub const TROPICAL_ZERO: f64 = 0.0;

/// The tropical semiring (max-plus algebra).
///
/// - Addition: a ⊕ b = max(a, b)
/// - Multiplication: a ⊗ b = a + b
/// - Additive identity: -∞
/// - Multiplicative identity: 0
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TropicalSemiring(pub f64);

impl TropicalSemiring {
    /// Create a tropical value. Use `None` or `NEG_INFINITY` for -∞.
    pub fn new(v: f64) -> Self {
        Self(v)
    }

    /// Tropical zero (-∞).
    pub fn zero() -> Self {
        Self(TROPICAL_NEG_INF)
    }

    /// Tropical one (0).
    pub fn one() -> Self {
        Self(TROPICAL_ZERO)
    }

    /// Is this the additive identity (-∞)?
    pub fn is_zero(&self) -> bool {
        self.0 == TROPICAL_NEG_INF || self.0.is_nan()
    }

    /// Tropical addition (max).
    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() { return *other; }
        if other.is_zero() { return *self; }
        Self(self.0.max(other.0))
    }

    /// Tropical multiplication (addition).
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        Self(self.0 + other.0)
    }

    /// Tropical power (scalar multiplication in conventional sense).
    pub fn pow(&self, n: usize) -> Self {
        if n == 0 { return Self::one(); }
        if self.is_zero() { return Self::zero(); }
        Self(self.0 * n as f64)
    }

    /// Tropical division (subtraction). Returns None if divisor is zero.
    pub fn div(&self, other: &Self) -> Option<Self> {
        if other.is_zero() { return None; }
        if self.is_zero() { return Some(Self::zero()); }
        Some(Self(self.0 - other.0))
    }
}

impl fmt::Display for TropicalSemiring {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() { write!(f, "-∞") } else { write!(f, "{}", self.0) }
    }
}

/// A matrix over the tropical semiring.
#[derive(Clone, Debug)]
pub struct TropicalMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<TropicalSemiring>,
}

impl TropicalMatrix {
    /// Create a tropical matrix filled with -∞.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![TropicalSemiring::zero(); rows * cols] }
    }

    /// Create a tropical identity matrix (0 on diagonal, -∞ elsewhere).
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n { m.set(i, i, TropicalSemiring::one()); }
        m
    }

    pub fn get(&self, i: usize, j: usize) -> TropicalSemiring {
        self.data[i * self.cols + j]
    }

    pub fn set(&mut self, i: usize, j: usize, v: TropicalSemiring) {
        self.data[i * self.cols + j] = v;
    }

    /// Tropical matrix multiplication: (A ⊗ B)_{ij} = max_k (A_{ik} + B_{kj}).
    pub fn matmul(&self, other: &Self) -> Self {
        assert_eq!(self.cols, other.rows);
        let mut result = Self::zeros(self.rows, other.cols);
        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut best = TropicalSemiring::zero();
                for k in 0..self.cols {
                    best = best.add(&self.get(i, k).mul(&other.get(k, j)));
                }
                result.set(i, j, best);
            }
        }
        result
    }

    /// Tropical matrix power.
    pub fn power(&self, n: usize) -> Self {
        if n == 0 { return Self::identity(self.rows); }
        if n == 1 { return self.clone(); }
        let mut result = self.clone();
        for _ in 1..n {
            result = result.matmul(self);
        }
        result
    }

    /// Kleene star: A* = I ⊕ A ⊕ A² ⊕ ... ⊕ Aⁿ (n = dimension).
    pub fn kleene_star(&self) -> Self {
        let n = self.rows;
        let mut result = Self::identity(n);
        let mut accum = self.clone();
        for _ in 0..n {
            result = tropical_matrix_add(&result, &accum);
            accum = accum.matmul(self);
        }
        result
    }
}

/// Element-wise tropical addition of two matrices.
pub fn tropical_matrix_add(a: &TropicalMatrix, b: &TropicalMatrix) -> TropicalMatrix {
    assert_eq!(a.rows, b.rows);
    assert_eq!(a.cols, b.cols);
    let mut result = TropicalMatrix::zeros(a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j).add(&b.get(i, j)));
        }
    }
    result
}

/// Shortest path computation using tropical matrix operations.
pub struct TropicalShortestPath;

impl TropicalShortestPath {
    /// Compute all-pairs longest paths via tropical matrix power (max-plus).
    /// For shortest paths, use min-plus instead.
    pub fn apsp(adj: &TropicalMatrix) -> TropicalMatrix {
        let n = adj.rows;
        let mut dist = adj.clone();
        for i in 0..n {
            dist.set(i, i, TropicalSemiring::one());
        }
        // Repeated squaring: dist = I ⊕ A ⊕ A² ⊕ ... ⊕ A^(n-1)
        let mut accum = dist.clone();
        let mut power = dist.clone();
        for _ in 0..n {
            power = power.matmul(&dist);
            for i in 0..n {
                for j in 0..n {
                    accum.set(i, j, accum.get(i, j).add(&power.get(i, j)));
                }
            }
        }
        accum
    }

    /// Compute all-pairs shortest paths using min-plus (min, +) semiring.
    pub fn shortest_paths(adj: &TropicalMatrix) -> TropicalMatrix {
        let n = adj.rows;
        // Floyd-Warshall style
        let mut dist = vec![vec![f64::INFINITY; n]; n];
        for i in 0..n { dist[i][i] = 0.0; }
        for i in 0..n {
            for j in 0..n {
                let v = adj.get(i, j);
                if !v.is_zero() {
                    dist[i][j] = v.0;
                }
            }
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    if dist[i][k] + dist[k][j] < dist[i][j] {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }
        let mut result = TropicalMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                if dist[i][j].is_finite() {
                    result.set(i, j, TropicalSemiring::new(dist[i][j]));
                }
            }
        }
        result
    }

    /// Single-source shortest paths from node `source`.
    pub fn sssp(adj: &TropicalMatrix, source: usize) -> Vec<TropicalSemiring> {
        let n = adj.rows;
        let mut dist = vec![TropicalSemiring::zero(); n];
        dist[source] = TropicalSemiring::one();

        // Tropical matrix-vector: repeated relaxation
        for _ in 0..n {
            let mut new_dist = dist.clone();
            for j in 0..n {
                for k in 0..n {
                    let candidate = dist[k].mul(&adj.get(k, j));
                    new_dist[j] = new_dist[j].add(&candidate);
                }
            }
            dist = new_dist;
        }
        dist
    }
}

/// Tropical eigenvalue (maximum cycle mean) computation.
pub struct TropicalEigenvalue;

impl TropicalEigenvalue {
    /// Compute the tropical eigenvalue (max cycle mean) via Karp's algorithm.
    ///
    /// λ = max_j min_k (Aⁿ[i][j] - Aᵏ[i][j]) / (n - k) for any fixed i.
    pub fn max_cycle_mean(adj: &TropicalMatrix) -> f64 {
        let n = adj.rows;
        if n == 0 { return f64::NEG_INFINITY; }

        // Compute A^0, A^1, ..., A^n from node 0
        let mut powers: Vec<Vec<TropicalSemiring>> = Vec::new();
        // A^0: identity row for node 0
        let mut row = vec![TropicalSemiring::zero(); n];
        row[0] = TropicalSemiring::one();
        powers.push(row);

        for p in 1..=n {
            let prev = &powers[p - 1];
            let mut new_row = vec![TropicalSemiring::zero(); n];
            for j in 0..n {
                for k in 0..n {
                    let candidate = prev[k].mul(&adj.get(k, j));
                    new_row[j] = new_row[j].add(&candidate);
                }
            }
            powers.push(new_row);
        }

        // Karp's formula: λ = max_j min_{k<n} ((Aⁿ)[j] - (Aᵏ)[j]) / (n - k)
        let mut lambda = f64::NEG_INFINITY;
        let an = &powers[n];
        for j in 0..n {
            if an[j].is_zero() { continue; }
            let mut min_ratio = f64::INFINITY;
            for k in 0..n {
                if powers[k][j].is_zero() { continue; }
                let diff = an[j].0 - powers[k][j].0;
                let ratio = diff / (n - k) as f64;
                if ratio < min_ratio {
                    min_ratio = ratio;
                }
            }
            if min_ratio != f64::INFINITY && min_ratio > lambda {
                lambda = min_ratio;
            }
        }
        lambda
    }
}

/// A tropical halfspace: { x : max(a₁ + x₁, a₂ + x₂, ..., aₙ + xₙ, b) achieved by b-side }.
/// Simplified: max_j(aⱼ + xⱼ) ≤ b  (or equivalently, the "≤" constraint in tropical sense).
#[derive(Clone, Debug)]
pub struct TropicalHalfspace {
    /// Coefficients a.
    pub coefficients: Vec<TropicalSemiring>,
    /// Bound b.
    pub bound: TropicalSemiring,
}

/// Tropical polyhedron: intersection of tropical halfspaces.
#[derive(Clone, Debug)]
pub struct TropicalPolyhedron {
    pub halfspaces: Vec<TropicalHalfspace>,
}

impl TropicalPolyhedron {
    /// Create an empty polyhedron (no constraints).
    pub fn new() -> Self {
        Self { halfspaces: Vec::new() }
    }

    /// Add a tropical halfspace constraint.
    pub fn add_halfspace(&mut self, coefficients: Vec<TropicalSemiring>, bound: TropicalSemiring) {
        self.halfspaces.push(TropicalHalfspace { coefficients, bound });
    }

    /// Check if a point satisfies all constraints.
    /// Constraint: max_j(aⱼ + xⱼ) ≤ b
    pub fn contains(&self, point: &[TropicalSemiring]) -> bool {
        for hs in &self.halfspaces {
            let mut max_val = TropicalSemiring::zero();
            for (j, a) in hs.coefficients.iter().enumerate() {
                if j < point.len() {
                    max_val = max_val.add(&a.mul(&point[j]));
                }
            }
            // Check: max_val ≤ bound
            if max_val.0 > hs.bound.0 && !max_val.is_zero() {
                return false;
            }
        }
        true
    }

    /// Number of constraints.
    pub fn num_constraints(&self) -> usize {
        self.halfspaces.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tropical_add() {
        let a = TropicalSemiring::new(3.0);
        let b = TropicalSemiring::new(5.0);
        assert_eq!(a.add(&b), TropicalSemiring::new(5.0));
        assert_eq!(b.add(&a), TropicalSemiring::new(5.0));
    }

    #[test]
    fn test_tropical_add_with_zero() {
        let a = TropicalSemiring::new(3.0);
        let z = TropicalSemiring::zero();
        assert_eq!(a.add(&z), a);
        assert_eq!(z.add(&a), a);
    }

    #[test]
    fn test_tropical_mul() {
        let a = TropicalSemiring::new(3.0);
        let b = TropicalSemiring::new(5.0);
        assert_eq!(a.mul(&b), TropicalSemiring::new(8.0));
    }

    #[test]
    fn test_tropical_mul_with_zero() {
        let a = TropicalSemiring::new(3.0);
        let z = TropicalSemiring::zero();
        assert!(a.mul(&z).is_zero());
        assert!(z.mul(&a).is_zero());
    }

    #[test]
    fn test_tropical_pow() {
        let a = TropicalSemiring::new(2.0);
        assert_eq!(a.pow(3), TropicalSemiring::new(6.0));
        assert_eq!(a.pow(0), TropicalSemiring::one());
    }

    #[test]
    fn test_tropical_div() {
        let a = TropicalSemiring::new(8.0);
        let b = TropicalSemiring::new(3.0);
        assert_eq!(a.div(&b), Some(TropicalSemiring::new(5.0)));
        assert!(TropicalSemiring::new(5.0).div(&TropicalSemiring::zero()).is_none());
    }

    #[test]
    fn test_tropical_identity_matrix() {
        let m = TropicalMatrix::identity(3);
        assert_eq!(m.get(0, 0), TropicalSemiring::one());
        assert!(m.get(0, 1).is_zero());
    }

    #[test]
    fn test_tropical_matrix_multiply() {
        let mut a = TropicalMatrix::zeros(2, 2);
        a.set(0, 0, TropicalSemiring::new(1.0));
        a.set(0, 1, TropicalSemiring::new(2.0));
        a.set(1, 0, TropicalSemiring::new(3.0));
        a.set(1, 1, TropicalSemiring::new(4.0));

        let mut b = TropicalMatrix::zeros(2, 2);
        b.set(0, 0, TropicalSemiring::new(5.0));
        b.set(0, 1, TropicalSemiring::new(6.0));
        b.set(1, 0, TropicalSemiring::new(7.0));
        b.set(1, 1, TropicalSemiring::new(8.0));

        let c = a.matmul(&b);
        // (0,0): max(1+5, 2+7) = max(6, 9) = 9
        assert_eq!(c.get(0, 0), TropicalSemiring::new(9.0));
        // (0,1): max(1+6, 2+8) = max(7, 10) = 10
        assert_eq!(c.get(0, 1), TropicalSemiring::new(10.0));
        // (1,0): max(3+5, 4+7) = max(8, 11) = 11
        assert_eq!(c.get(1, 0), TropicalSemiring::new(11.0));
    }

    #[test]
    fn test_tropical_matrix_power() {
        let mut a = TropicalMatrix::zeros(2, 2);
        a.set(0, 0, TropicalSemiring::new(1.0));
        a.set(0, 1, TropicalSemiring::new(2.0));
        a.set(1, 0, TropicalSemiring::new(3.0));
        a.set(1, 1, TropicalSemiring::new(1.0));

        let a2 = a.power(2);
        // A²[0][0] = max(1+1, 2+3) = max(2, 5) = 5
        assert_eq!(a2.get(0, 0), TropicalSemiring::new(5.0));
    }

    #[test]
    fn test_shortest_path_triangle() {
        let mut adj = TropicalMatrix::zeros(3, 3);
        adj.set(0, 1, TropicalSemiring::new(2.0));
        adj.set(1, 2, TropicalSemiring::new(3.0));
        adj.set(0, 2, TropicalSemiring::new(10.0));

        let dist = TropicalShortestPath::shortest_paths(&adj);
        // Shortest 0→2 should be min(10, 2+3) = 5
        assert_eq!(dist.get(0, 2), TropicalSemiring::new(5.0));
        // Self distances should be 0
        assert_eq!(dist.get(0, 0), TropicalSemiring::new(0.0));
    }

    #[test]
    fn test_sssp() {
        let mut adj = TropicalMatrix::zeros(3, 3);
        adj.set(0, 1, TropicalSemiring::new(2.0));
        adj.set(1, 2, TropicalSemiring::new(3.0));

        let dist = TropicalShortestPath::sssp(&adj, 0);
        assert_eq!(dist[0], TropicalSemiring::new(0.0));
        assert_eq!(dist[1], TropicalSemiring::new(2.0));
        assert_eq!(dist[2], TropicalSemiring::new(5.0));
    }

    #[test]
    fn test_max_cycle_mean() {
        let mut adj = TropicalMatrix::zeros(2, 2);
        adj.set(0, 0, TropicalSemiring::new(3.0));
        adj.set(0, 1, TropicalSemiring::new(1.0));
        adj.set(1, 0, TropicalSemiring::new(2.0));
        adj.set(1, 1, TropicalSemiring::new(5.0));

        let lambda = TropicalEigenvalue::max_cycle_mean(&adj);
        // Max cycle mean = max(3, 5, (1+2)/2) = max(3, 5, 1.5) = 5
        assert!((lambda - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_tropical_polyhedron_contains() {
        let mut poly = TropicalPolyhedron::new();
        // max(x₁, x₂) ≤ 5
        poly.add_halfspace(
            vec![TropicalSemiring::new(0.0), TropicalSemiring::new(0.0)],
            TropicalSemiring::new(5.0),
        );
        assert!(poly.contains(&[TropicalSemiring::new(3.0), TropicalSemiring::new(4.0)]));
        assert!(!poly.contains(&[TropicalSemiring::new(6.0), TropicalSemiring::new(3.0)]));
    }

    #[test]
    fn test_tropical_polyhedron_empty() {
        let poly = TropicalPolyhedron::new();
        assert!(poly.contains(&[TropicalSemiring::new(100.0)]));
        assert_eq!(poly.num_constraints(), 0);
    }
}
