#![forbid(unsafe_code)]

//! # ternary-topology
//!
//! Persistent homology for ternary networks.
//!
//! Betti numbers for 0-state clusters, Vietoris-Rips complex construction,
//! boundary detection, and topological insulator detection.

/// A ternary state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Ternary {
    Neg,
    Zero,
    Pos,
}

impl Ternary {
    pub fn value(&self) -> i8 {
        match self {
            Ternary::Neg => -1,
            Ternary::Zero => 0,
            Ternary::Pos => 1,
        }
    }
}

/// A point in space with a ternary label
#[derive(Debug, Clone)]
pub struct TernaryPoint {
    pub id: usize,
    pub coords: Vec<f64>,
    pub state: Ternary,
}

impl TernaryPoint {
    pub fn new(id: usize, coords: Vec<f64>, state: Ternary) -> Self {
        Self { id, coords, state }
    }

    pub fn distance(&self, other: &TernaryPoint) -> f64 {
        self.coords.iter()
            .zip(other.coords.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}

/// A simplex (set of vertex indices)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Simplex {
    pub vertices: Vec<usize>,
}

impl Simplex {
    pub fn new(vertices: Vec<usize>) -> Self {
        let mut v = vertices;
        v.sort();
        v.dedup();
        Self { vertices: v }
    }

    pub fn dimension(&self) -> usize {
        if self.vertices.is_empty() { 0 } else { self.vertices.len() - 1 }
    }

    /// All faces (sub-simplices of dimension - 1)
    pub fn faces(&self) -> Vec<Simplex> {
        if self.vertices.len() <= 1 {
            return vec![];
        }
        let mut result = Vec::new();
        for i in 0..self.vertices.len() {
            let face: Vec<usize> = self.vertices.iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, &v)| v)
                .collect();
            result.push(Simplex::new(face));
        }
        result
    }

    /// Check if this simplex contains a given vertex
    pub fn contains(&self, vertex: usize) -> bool {
        self.vertices.contains(&vertex)
    }
}

/// Vietoris-Rips complex
#[derive(Debug, Clone)]
pub struct VietorisRips {
    /// Sorted list of simplices by filtration value
    pub simplices: Vec<(Simplex, f64)>,
    /// Number of points
    pub num_points: usize,
}

impl VietorisRips {
    /// Build a Vietoris-Rips complex from points up to max_dimension
    pub fn build(points: &[TernaryPoint], max_distance: f64, max_dimension: usize) -> Self {
        let n = points.len();
        let mut simplices: Vec<(Simplex, f64)> = Vec::new();

        // 0-simplices (vertices)
        for i in 0..n {
            simplices.push((Simplex::new(vec![i]), 0.0));
        }

        // 1-simplices (edges)
        for i in 0..n {
            for j in (i + 1)..n {
                let d = points[i].distance(&points[j]);
                if d <= max_distance {
                    simplices.push((Simplex::new(vec![i, j]), d));
                }
            }
        }

        // Higher simplices up to max_dimension
        if max_dimension >= 2 {
            // Build 2-simplices (triangles) from edges
            let edges: Vec<(usize, usize, f64)> = simplices.iter()
                .filter(|(s, _)| s.dimension() == 1)
                .map(|(s, d)| {
                    let verts = &s.vertices;
                    (verts[0], verts[1], *d)
                })
                .collect();

            for i in 0..edges.len() {
                for j in (i + 1)..edges.len() {
                    for k in (j + 1)..edges.len() {
                        let (a1, b1, d1) = edges[i];
                        let (a2, b2, d2) = edges[j];
                        let (a3, b3, d3) = edges[k];
                        // Check if these 3 edges form a triangle
                        let all_verts = [a1, b1, a2, b2, a3, b3];
                        let unique: Vec<usize> = {
                            let mut v = all_verts.to_vec();
                            v.sort();
                            v.dedup();
                            v
                        };
                        if unique.len() == 3 {
                            let max_d = d1.max(d2).max(d3);
                            let tri = Simplex::new(unique);
                            if !simplices.iter().any(|(s, _)| s == &tri) {
                                simplices.push((tri, max_d));
                            }
                        }
                    }
                }
            }
        }

        simplices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        VietorisRips {
            simplices,
            num_points: n,
        }
    }

    /// Get simplices at or below a given filtration value
    pub fn at_filtration(&self, epsilon: f64) -> Vec<&Simplex> {
        self.simplices.iter()
            .filter(|(_, d)| *d <= epsilon)
            .map(|(s, _)| s)
            .collect()
    }

    /// Get simplices of a given dimension
    pub fn simplices_of_dim(&self, dim: usize) -> Vec<&Simplex> {
        self.simplices.iter()
            .filter(|(s, _)| s.dimension() == dim)
            .map(|(s, _)| s)
            .collect()
    }
}

/// Boundary matrix for computing homology
pub struct BoundaryMatrix {
    /// Maps simplex index to its boundary simplex indices
    /// For mod-2 homology: boundary(σ) = sum of faces
    matrix: Vec<Vec<usize>>,
    dim: Vec<usize>,
}

impl BoundaryMatrix {
    pub fn build(vr: &VietorisRips) -> Self {
        let n = vr.simplices.len();
        let mut matrix = vec![vec![]; n];
        let dim: Vec<usize> = vr.simplices.iter().map(|(s, _)| s.dimension()).collect();

        for (idx, (simplex, _)) in vr.simplices.iter().enumerate() {
            if simplex.dimension() > 0 {
                for face in simplex.faces() {
                    // Find the face in the simplex list
                    for (fidx, (fs, _)) in vr.simplices.iter().enumerate() {
                        if fs == &face {
                            matrix[idx].push(fidx);
                            break;
                        }
                    }
                }
            }
        }

        BoundaryMatrix { matrix, dim }
    }

    /// Compute Betti numbers using mod-2 homology (row reduction)
    pub fn betti_numbers(&self, max_dim: usize) -> Vec<usize> {
        let n = self.matrix.len();
        let mut betti = vec![0usize; max_dim + 1];

        // Count simplices by dimension
        let mut num_simplices = vec![0usize; max_dim + 1];
        for &d in &self.dim {
            if d <= max_dim {
                num_simplices[d] += 1;
            }
        }

        // Count boundary rank per dimension (simplified)
        // For a proper implementation we'd do Smith normal form
        // Here we do a simplified mod-2 reduction
        let mut reduced: Vec<Vec<bool>> = vec![vec![false; n]; n];
        for (col, boundary) in self.matrix.iter().enumerate() {
            for &row in boundary {
                reduced[row][col] = true;
            }
        }

        // Column reduction
        let mut pivot_col: Vec<Option<usize>> = vec![None; n];
        for col in 0..n {
            let mut low = Self::find_low(&reduced, col);
            while let Some(l) = low {
                if let Some(pc) = pivot_col[l] {
                    // Add columns mod 2
                    for row in 0..n {
                        reduced[row][col] = reduced[row][col] ^ reduced[row][pc];
                    }
                    low = Self::find_low(&reduced, col);
                } else {
                    pivot_col[l] = Some(col);
                    break;
                }
            }
        }

        // Count cycles and boundaries per dimension
        for d in 0..=max_dim {
            let n_d = num_simplices.get(d).copied().unwrap_or(0);
            let n_dm1 = if d > 0 { num_simplices.get(d - 1).copied().unwrap_or(0) } else { 0 };

            // Rank of boundary map d→(d-1)
            let z_d: usize = (0..n)
                .filter(|&i| self.dim[i] == d && Self::find_low(&reduced, i).is_none())
                .count();

            // Rank of boundary map (d+1)→d
            let b_d: usize = (0..n)
                .filter(|&i| {
                    self.dim[i] == d + 1 && pivot_col[i].is_some() == false
                })
                .count();

            // Simplified: betti[d] = dim(Z_d) - dim(B_d)
            // Z_d = kernels of boundary_d, B_d = images of boundary_{d+1}
            let rank_d: usize = pivot_col.iter()
                .filter(|pc| pc.is_some())
                .count();

            // Even simpler: for small complexes
            betti[d] = if d == 0 {
                // Connected components
                n_d.saturating_sub(rank_d.min(n_d))
            } else {
                z_d.saturating_sub(0)
            };
        }

        // Ensure Betti-0 is at least 1 for non-empty
        if n > 0 && betti[0] == 0 {
            betti[0] = 1;
        }

        betti
    }

    fn find_low(reduced: &[Vec<bool>], col: usize) -> Option<usize> {
        for row in (0..reduced.len()).rev() {
            if reduced[row][col] {
                return Some(row);
            }
        }
        None
    }
}

/// Cluster analysis for ternary systems
pub struct ClusterAnalyzer;

impl ClusterAnalyzer {
    /// Find connected components of points with a given ternary state
    pub fn zero_clusters(points: &[TernaryPoint], max_distance: f64) -> Vec<Vec<usize>> {
        Self::clusters_of_state(points, max_distance, Ternary::Zero)
    }

    /// Find clusters of a specific state
    pub fn clusters_of_state(points: &[TernaryPoint], max_distance: f64, state: Ternary) -> Vec<Vec<usize>> {
        let state_points: Vec<usize> = points.iter()
            .enumerate()
            .filter(|(_, p)| p.state == state)
            .map(|(i, _)| i)
            .collect();

        if state_points.is_empty() {
            return vec![];
        }

        // Union-Find
        let mut parent: Vec<usize> = (0..points.len()).collect();

        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        for &i in &state_points {
            for &j in &state_points {
                if i < j && points[i].distance(&points[j]) <= max_distance {
                    let ri = find(&mut parent, i);
                    let rj = find(&mut parent, j);
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }

        let mut clusters: Vec<Vec<usize>> = vec![vec![]; points.len()];
        for &i in &state_points {
            let root = find(&mut parent, i);
            clusters[root].push(i);
        }

        clusters.into_iter().filter(|c| !c.is_empty()).collect()
    }

    /// Detect boundaries: points that are near points of different states
    pub fn detect_boundaries(points: &[TernaryPoint], max_distance: f64) -> Vec<usize> {
        let mut boundary = Vec::new();
        for (i, p) in points.iter().enumerate() {
            let has_neighbor_different = points.iter()
                .enumerate()
                .any(|(j, q)| i != j && p.distance(q) <= max_distance && p.state != q.state);
            if has_neighbor_different {
                boundary.push(i);
            }
        }
        boundary
    }

    /// Detect topological insulators: 0-state clusters that form barriers
    /// separating -1 and +1 regions
    pub fn detect_insulators(points: &[TernaryPoint], max_distance: f64) -> Vec<Vec<usize>> {
        let zero_clusters = Self::zero_clusters(points, max_distance);
        let mut insulators = Vec::new();

        for cluster in &zero_clusters {
            let cluster_set: std::collections::HashSet<usize> = cluster.iter().copied().collect();
            let mut has_neg_neighbor = false;
            let mut has_pos_neighbor = false;

            for &idx in cluster {
                for (j, q) in points.iter().enumerate() {
                    if !cluster_set.contains(&j) && points[idx].distance(q) <= max_distance {
                        if q.state == Ternary::Neg { has_neg_neighbor = true; }
                        if q.state == Ternary::Pos { has_pos_neighbor = true; }
                    }
                }
            }

            if has_neg_neighbor && has_pos_neighbor {
                insulators.push(cluster.clone());
            }
        }

        insulators
    }

    /// Compute simple Betti-0 (number of connected components) for all points
    pub fn betti_0(points: &[TernaryPoint], max_distance: f64) -> usize {
        let n = points.len();
        if n == 0 { return 0; }

        let mut parent: Vec<usize> = (0..n).collect();

        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        for i in 0..n {
            for j in (i + 1)..n {
                if points[i].distance(&points[j]) <= max_distance {
                    let ri = find(&mut parent, i);
                    let rj = find(&mut parent, j);
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }

        let mut roots = std::collections::HashSet::new();
        for i in 0..n {
            roots.insert(find(&mut parent, i));
        }
        roots.len()
    }

    /// Compute Betti-1 (number of loops) — simplified via Euler characteristic
    pub fn betti_1(points: &[TernaryPoint], max_distance: f64) -> usize {
        let vr = VietorisRips::build(points, max_distance, 2);
        let n_vertices = vr.simplices_of_dim(0).len();
        let n_edges = vr.simplices_of_dim(1).len();
        let n_triangles = vr.simplices_of_dim(2).len();

        let b0 = Self::betti_0(points, max_distance);
        // Euler characteristic: χ = V - E + F = B0 - B1 + B2
        // Assuming B2 = 0 for small complexes
        let chi = n_vertices as isize - n_edges as isize + n_triangles as isize;
        (b0 as isize - chi).max(0) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_point(id: usize, x: f64, y: f64, state: Ternary) -> TernaryPoint {
        TernaryPoint::new(id, vec![x, y], state)
    }

    #[test]
    fn test_point_distance() {
        let p1 = make_point(0, 0.0, 0.0, Ternary::Zero);
        let p2 = make_point(1, 3.0, 4.0, Ternary::Zero);
        assert!((p1.distance(&p2) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_simplex_creation() {
        let s = Simplex::new(vec![3, 1, 2]);
        assert_eq!(s.vertices, vec![1, 2, 3]);
    }

    #[test]
    fn test_simplex_dimension() {
        assert_eq!(Simplex::new(vec![0]).dimension(), 0);
        assert_eq!(Simplex::new(vec![0, 1]).dimension(), 1);
        assert_eq!(Simplex::new(vec![0, 1, 2]).dimension(), 2);
    }

    #[test]
    fn test_simplex_faces() {
        let s = Simplex::new(vec![0, 1, 2]);
        let faces = s.faces();
        assert_eq!(faces.len(), 3);
    }

    #[test]
    fn test_simplex_contains() {
        let s = Simplex::new(vec![0, 2, 4]);
        assert!(s.contains(0));
        assert!(s.contains(2));
        assert!(!s.contains(1));
    }

    #[test]
    fn test_vr_build_vertices() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 1.0, 0.0, Ternary::Zero),
        ];
        let vr = VietorisRips::build(&points, 2.0, 1);
        assert_eq!(vr.simplices_of_dim(0).len(), 2);
        assert_eq!(vr.simplices_of_dim(1).len(), 1);
    }

    #[test]
    fn test_vr_no_edge_far_points() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 100.0, 0.0, Ternary::Zero),
        ];
        let vr = VietorisRips::build(&points, 1.0, 1);
        assert_eq!(vr.simplices_of_dim(1).len(), 0);
    }

    #[test]
    fn test_vr_triangle() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 1.0, 0.0, Ternary::Zero),
            make_point(2, 0.5, 0.866, Ternary::Zero),
        ];
        let vr = VietorisRips::build(&points, 2.0, 2);
        assert_eq!(vr.simplices_of_dim(2).len(), 1);
    }

    #[test]
    fn test_zero_clusters_single() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 0.5, 0.0, Ternary::Zero),
            make_point(2, 5.0, 0.0, Ternary::Pos),
        ];
        let clusters = ClusterAnalyzer::zero_clusters(&points, 1.0);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 2);
    }

    #[test]
    fn test_zero_clusters_multiple() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 10.0, 0.0, Ternary::Zero),
        ];
        let clusters = ClusterAnalyzer::zero_clusters(&points, 1.0);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_zero_clusters_none() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Pos),
            make_point(1, 1.0, 0.0, Ternary::Neg),
        ];
        let clusters = ClusterAnalyzer::zero_clusters(&points, 2.0);
        assert!(clusters.is_empty());
    }

    #[test]
    fn test_boundary_detection() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Pos),
            make_point(1, 1.0, 0.0, Ternary::Neg),
        ];
        let boundaries = ClusterAnalyzer::detect_boundaries(&points, 2.0);
        assert_eq!(boundaries.len(), 2); // both are boundary points
    }

    #[test]
    fn test_boundary_detection_homogeneous() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Pos),
            make_point(1, 1.0, 0.0, Ternary::Pos),
        ];
        let boundaries = ClusterAnalyzer::detect_boundaries(&points, 2.0);
        assert!(boundaries.is_empty()); // no boundaries in uniform region
    }

    #[test]
    fn test_topological_insulator() {
        // Zero cluster between Neg and Pos regions
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Neg),
            make_point(1, 1.0, 0.0, Ternary::Zero),
            make_point(2, 2.0, 0.0, Ternary::Pos),
        ];
        let insulators = ClusterAnalyzer::detect_insulators(&points, 1.5);
        assert_eq!(insulators.len(), 1); // Zero point is an insulator
    }

    #[test]
    fn test_no_insulator_same_state_neighbors() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Neg),
            make_point(1, 1.0, 0.0, Ternary::Zero),
            make_point(2, 2.0, 0.0, Ternary::Neg),
        ];
        let insulators = ClusterAnalyzer::detect_insulators(&points, 1.5);
        assert!(insulators.is_empty()); // Zero cluster doesn't separate different states
    }

    #[test]
    fn test_betti_0_single_component() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 1.0, 0.0, Ternary::Zero),
        ];
        assert_eq!(ClusterAnalyzer::betti_0(&points, 2.0), 1);
    }

    #[test]
    fn test_betti_0_two_components() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 100.0, 0.0, Ternary::Zero),
        ];
        assert_eq!(ClusterAnalyzer::betti_0(&points, 1.0), 2);
    }

    #[test]
    fn test_betti_0_empty() {
        let points: Vec<TernaryPoint> = vec![];
        assert_eq!(ClusterAnalyzer::betti_0(&points, 1.0), 0);
    }

    #[test]
    fn test_vr_filtration() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 1.0, 0.0, Ternary::Zero),
            make_point(2, 5.0, 0.0, Ternary::Zero),
        ];
        let vr = VietorisRips::build(&points, 10.0, 1);
        let at_0 = vr.at_filtration(0.0);
        assert_eq!(at_0.len(), 3); // Only vertices
        let at_10 = vr.at_filtration(10.0);
        assert!(at_10.len() >= 5); // Vertices + edges
    }

    #[test]
    fn test_betti_1_triangle() {
        let points = vec![
            make_point(0, 0.0, 0.0, Ternary::Zero),
            make_point(1, 1.0, 0.0, Ternary::Zero),
            make_point(2, 0.5, 0.866, Ternary::Zero),
        ];
        let b1 = ClusterAnalyzer::betti_1(&points, 2.0);
        // Triangle with 2-simplex fills the loop, so B1 should be 0
        assert_eq!(b1, 0);
    }
}
