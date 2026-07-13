//! # ternary-topology
//!
//! Self-organizing network topology via ternary adjacency.
//! Nodes form edges based on ternary affinity {-1, 0, +1}.

use std::collections::HashMap;

/// A node in the topology.
#[derive(Debug, Clone)]
pub struct Node {
    pub id: u32,
    pub ternary_state: Vec<i8>,
    pub neighbors: Vec<u32>,
}

impl Node {
    pub fn new(id: u32, dim: usize) -> Self {
        Self {
            id,
            ternary_state: vec![0; dim],
            neighbors: Vec::new(),
        }
    }
}

/// Ternary affinity between two nodes based on their state vectors.
pub fn affinity(a: &[i8], b: &[i8]) -> i8 {
    let sum: i32 = a.iter().zip(b).map(|(&x, &y)| x as i32 * y as i32).sum();
    if sum > 0 {
        1
    } else if sum < 0 {
        -1
    } else {
        0
    }
}

/// A self-organizing topology.
pub struct TernaryTopology {
    pub nodes: HashMap<u32, Node>,
    pub edges: Vec<(u32, u32, i8)>, // (from, to, affinity)
    pub dimension: usize,
}

impl TernaryTopology {
    pub fn new(dimension: usize) -> Self {
        Self {
            nodes: HashMap::new(),
            edges: Vec::new(),
            dimension,
        }
    }

    /// Add a node with a ternary state vector.
    ///
    /// Panics if `state.len()` does not equal the topology's `dimension`.
    /// The fixed-dimension invariant is what guarantees that `affinity` and
    /// `evolve` are well-defined for every pair of nodes; accepting a
    /// mismatched length would otherwise be truncated silently by `affinity`
    /// and could index out of bounds inside `evolve`.
    pub fn add_node(&mut self, id: u32, state: Vec<i8>) {
        assert_eq!(
            state.len(),
            self.dimension,
            "node {} state length {} does not match topology dimension {}",
            id,
            state.len(),
            self.dimension,
        );
        self.nodes.insert(
            id,
            Node {
                id,
                ternary_state: state,
                neighbors: Vec::new(),
            },
        );
    }

    /// Compute all edges based on current states.
    ///
    /// Edges are emitted in ascending order of node id, with each undirected
    /// edge normalized so the lower id comes first. Iterating the node set in
    /// sorted order makes the resulting edge list (and therefore
    /// `positive_edges` / `negative_edges`) deterministic and reproducible.
    pub fn compute_edges(&mut self) {
        self.edges.clear();
        let mut ids: Vec<u32> = self.nodes.keys().cloned().collect();
        ids.sort_unstable();
        for i in 0..ids.len() {
            for j in (i + 1)..ids.len() {
                let a = &self.nodes[&ids[i]];
                let b = &self.nodes[&ids[j]];
                let aff = affinity(&a.ternary_state, &b.ternary_state);
                if aff != 0 {
                    self.edges.push((ids[i], ids[j], aff));
                }
            }
        }
        // Update neighbor lists
        for node in self.nodes.values_mut() {
            node.neighbors.clear();
        }
        for &(a, b, aff) in &self.edges {
            if aff > 0 {
                self.nodes.get_mut(&a).unwrap().neighbors.push(b);
                self.nodes.get_mut(&b).unwrap().neighbors.push(a);
            }
        }
    }

    /// Evolve: each node moves toward the average of positive-affinity neighbors.
    pub fn evolve(&mut self, steps: usize) -> Vec<usize> {
        let mut edge_counts = Vec::with_capacity(steps);
        for _ in 0..steps {
            self.compute_edges();
            let edge_count = self.edges.len();
            edge_counts.push(edge_count);

            // Compute new states
            let mut new_states: HashMap<u32, Vec<i8>> = HashMap::new();
            for (&id, node) in &self.nodes {
                let mut avg = vec![0f64; self.dimension];
                let mut count = 0;
                for &nid in &node.neighbors {
                    if let Some(n) = self.nodes.get(&nid) {
                        for (i, &v) in n.ternary_state.iter().enumerate() {
                            avg[i] += v as f64;
                        }
                        count += 1;
                    }
                }
                if count > 0 {
                    let new_state: Vec<i8> = avg
                        .iter()
                        .map(|&v| {
                            let mean = v / count as f64;
                            if mean > 0.3 {
                                1
                            } else if mean < -0.3 {
                                -1
                            } else {
                                0
                            }
                        })
                        .collect();
                    new_states.insert(id, new_state);
                }
            }
            for (id, state) in new_states {
                self.nodes.get_mut(&id).unwrap().ternary_state = state;
            }
        }
        edge_counts
    }

    pub fn positive_edges(&self) -> Vec<(u32, u32)> {
        self.edges
            .iter()
            .filter(|(_, _, a)| *a > 0)
            .map(|&(a, b, _)| (a, b))
            .collect()
    }

    pub fn negative_edges(&self) -> Vec<(u32, u32)> {
        self.edges
            .iter()
            .filter(|(_, _, a)| *a < 0)
            .map(|&(a, b, _)| (a, b))
            .collect()
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_affinity_positive() {
        assert_eq!(affinity(&[1, 1, 1], &[1, 1, 1]), 1);
    }

    #[test]
    fn test_affinity_negative() {
        assert_eq!(affinity(&[1, 1, 1], &[-1, -1, -1]), -1);
    }

    #[test]
    fn test_affinity_neutral() {
        assert_eq!(affinity(&[0, 0, 0], &[1, -1, 1]), 0);
        assert_eq!(affinity(&[1, -1, 0], &[0, 0, 1]), 0);
    }

    #[test]
    fn test_compute_edges() {
        let mut topo = TernaryTopology::new(3);
        topo.add_node(0, vec![1, 1, 1]);
        topo.add_node(1, vec![1, 1, 1]);
        topo.add_node(2, vec![-1, -1, -1]);
        topo.compute_edges();
        assert_eq!(topo.positive_edges().len(), 1); // 0↔1
        assert_eq!(topo.negative_edges().len(), 2); // 0↔2, 1↔2
    }

    #[test]
    fn test_edges_are_deterministic_and_normalized() {
        // The undirected edge (u,v) is always reported with u < v, and the
        // full edge list is in ascending order regardless of insertion order.
        let mut topo = TernaryTopology::new(3);
        topo.add_node(2, vec![-1, -1, -1]);
        topo.add_node(0, vec![1, 1, 1]);
        topo.add_node(1, vec![1, 1, 1]);
        topo.compute_edges();
        assert_eq!(topo.positive_edges(), vec![(0, 1)]);
        assert_eq!(topo.negative_edges(), vec![(0, 2), (1, 2)]);
    }

    #[test]
    fn test_quickstart_edge_correctness() {
        // Mirrors the README Quick Start. Verifies the actual edge set so the
        // documented example can never silently drift again.
        let mut topo = TernaryTopology::new(4);
        topo.add_node(0, vec![1, 1, 1, 1]);
        topo.add_node(1, vec![1, 1, 0, 1]);
        topo.add_node(2, vec![-1, -1, -1, 0]);
        topo.add_node(3, vec![-1, -1, 0, -1]);
        topo.compute_edges();
        assert_eq!(topo.positive_edges(), vec![(0, 1), (2, 3)]);
        assert_eq!(topo.negative_edges(), vec![(0, 2), (0, 3), (1, 2), (1, 3)]);
    }

    #[test]
    fn test_evolution_converges() {
        let mut topo = TernaryTopology::new(4);
        for i in 0..6 {
            topo.add_node(i, vec![if i < 3 { 1 } else { -1 }; 4]);
        }
        let counts = topo.evolve(5);
        assert_eq!(counts.len(), 5);
        // Two clusters of three aligned nodes: 3 positive edges inside each
        // cluster and 3*3 = 9 negative edges between them. The configuration
        // is a fixed point of the evolution rule, so the count is constant.
        assert!(
            counts.iter().all(|&c| c == 15),
            "edge counts were {counts:?}"
        );
        // The evolution rule pulls each node toward its identical neighbors,
        // so the two clusters stay at their respective poles.
        for i in 0..3 {
            assert_eq!(topo.nodes[&i].ternary_state, vec![1, 1, 1, 1]);
        }
        for i in 3..6 {
            assert_eq!(topo.nodes[&i].ternary_state, vec![-1, -1, -1, -1]);
        }
    }

    #[test]
    fn test_neutral_node_is_isolated() {
        // A zero state is orthogonal to everyone, so it forms no edges and is
        // never pulled during evolution (the README "boundary node" case).
        let mut topo = TernaryTopology::new(3);
        topo.add_node(0, vec![1, 1, 1]);
        topo.add_node(1, vec![1, 1, 1]);
        topo.add_node(2, vec![0, 0, 0]);
        topo.compute_edges();
        assert_eq!(topo.positive_edges(), vec![(0, 1)]);
        assert!(topo.negative_edges().is_empty());
        assert!(topo.nodes[&2].neighbors.is_empty());
        topo.evolve(3);
        assert_eq!(topo.nodes[&2].ternary_state, vec![0, 0, 0]);
    }

    #[test]
    #[should_panic(expected = "does not match topology dimension")]
    fn test_dimension_mismatch_panics() {
        // The "fixed dimension" invariant must be enforced: a node whose
        // state length disagrees with the topology dimension used to be
        // accepted silently, causing affinity to truncate and evolve to risk
        // an out-of-bounds index. It must now panic.
        let mut topo = TernaryTopology::new(4);
        topo.add_node(0, vec![1, 1, 1, 1]);
        topo.add_node(1, vec![1, 1]); // wrong length
    }

    #[test]
    fn test_neighbors_connected() {
        let mut topo = TernaryTopology::new(3);
        topo.add_node(0, vec![1, 1, 1]);
        topo.add_node(1, vec![1, 1, 1]);
        topo.compute_edges();
        assert!(topo.nodes[&0].neighbors.contains(&1));
        assert!(topo.nodes[&1].neighbors.contains(&0));
    }

    #[test]
    fn test_empty_topology() {
        let topo = TernaryTopology::new(4);
        assert_eq!(topo.node_count(), 0);
        assert_eq!(topo.positive_edges().len(), 0);
    }
}
