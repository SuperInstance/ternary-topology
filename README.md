# ternary-topology

**Self-organizing networks where nodes find their neighbors through ternary affinity. State vectors of {-1, 0, +1} drive edge formation, evolution, and convergence — no central coordinator needed.**

## Why This Exists

Most network topologies are designed: a human or algorithm decides who connects to whom. But in biological neural networks, social networks, and distributed systems, connections form organically — nodes discover each other based on local affinity, not global planning.

This crate implements self-organizing topology formation using ternary state vectors. Each node carries a state vector of {-1, 0, +1} values. The *affinity* between two nodes — whether they attract (+1), repel (-1), or ignore (0) each other — is computed from their state vectors. Edges form where affinity is nonzero. Over time, nodes evolve their states toward their positive-affinity neighbors, and the topology converges to a stable structure.

No coordinator. No global state. Just local affinity driving global organization.

## The Key Insight

The affinity function is the dot product projected to ternary:

```rust
pub fn affinity(a: &[i8], b: &[i8]) -> i8 {
    let sum: i32 = a.iter().zip(b).map(|(&x, &y)| x as i32 * y as i32).sum();
    if sum > 0 { 1 } else if sum < 0 { -1 } else { 0 }
}
```

This is ternary correlation: the sign of the inner product tells you whether two state vectors are aligned (+1), opposed (-1), or orthogonal (0). Aligned nodes form positive edges and pull each other closer. Opposed nodes form negative edges and push apart. Orthogonal nodes don't interact at all — they're invisible to each other.

What makes this powerful is the *evolution rule*: after computing edges, each node moves its state toward the mean of its positive-affinity neighbors (thresholded back to {-1, 0, +1}). This creates a feedback loop:

1. State determines affinity → affinity determines edges → edges determine neighbors
2. Neighbors determine evolution → evolution changes state → repeat

The system converges when nodes have settled into clusters of aligned states, with negative edges forming boundaries between clusters.

## Quick Start

```rust
use ternary_topology::*;

// Create a topology with 4-dimensional ternary state vectors
let mut topo = TernaryTopology::new(4);

// Add nodes with initial states
topo.add_node(0, vec![1, 1, 1, 1]);   // Cluster A
topo.add_node(1, vec![1, 1, 0, 1]);   // Cluster A (similar)
topo.add_node(2, vec![-1, -1, -1, 0]); // Cluster B (opposed)
topo.add_node(3, vec![-1, -1, 0, -1]); // Cluster B (similar)

// Compute edges based on affinity
topo.compute_edges();
let positive = topo.positive_edges();  // Edges with affinity +1
let negative = topo.negative_edges();  // Edges with affinity -1
println!("Positive connections: {:?}", positive);  // 0↔1, 2↔3
println!("Negative connections: {:?}", negative);  // 0↔2, 0↔3, 1↔2, 1↔3

// Evolve: nodes move toward their positive neighbors
let edge_counts = topo.evolve(5);
// edge_counts tracks how many edges exist at each step
// As nodes converge, edge counts may change
```

## Architecture

### Core Types

| Type | Description |
|------|-------------|
| `Node` | A node with ID, ternary state vector, and neighbor list |
| `TernaryTopology` | The network: set of nodes, computed edges, and dimension |

### Affinity Function

```
affinity(a, b) = sign(a · b)

Examples:
  affinity([1,1,1], [1,1,1]) = +1  (perfectly aligned)
  affinity([1,1,1], [-1,-1,-1]) = -1  (perfectly opposed)
  affinity([1,-1,0], [0,0,1]) = 0  (orthogonal)
  affinity([1,-1,1], [1,-1,-1]) = +1  (2/3 aligned → positive sum)
```

### TernaryTopology Methods

| Method | Description |
|--------|-------------|
| `new(dimension)` | Create empty topology with given state dimension |
| `add_node(id, state)` | Add node with ternary state vector |
| `compute_edges()` | Compute all pairwise affinities, build edge list (ascending, deterministic) |
| `evolve(steps)` | Run N evolution steps, return edge counts per step |
| `positive_edges()` | Edges with affinity +1 (attraction) |
| `negative_edges()` | Edges with affinity -1 (repulsion) |
| `node_count()` | Number of nodes |

### Evolution Mechanics

Each evolution step:

1. **Recompute edges** — affinity recalculated from current states
2. **Gather neighbor states** — for each node, collect states of positive-affinity neighbors
3. **Compute mean** — element-wise average of neighbor states
4. **Ternarize** — threshold each dimension: `> 0.3 → +1`, `< -0.3 → -1`, else `0`
5. **Update state** — node adopts the ternarized mean

The 0.3 threshold prevents small fluctuations from creating state changes. A dimension only flips when there's a clear signal from the majority of neighbors.

## Real-World Example: Cluster Formation

```rust
use ternary_topology::*;

let mut topo = TernaryTopology::new(3);

// Two clusters of aligned nodes, plus a boundary node
// Cluster A: nodes 0-2
topo.add_node(0, vec![1, 1, 1]);
topo.add_node(1, vec![1, 1, 0]);
topo.add_node(2, vec![1, 0, 1]);

// Cluster B: nodes 3-5
topo.add_node(3, vec![-1, -1, -1]);
topo.add_node(4, vec![-1, -1, 0]);
topo.add_node(5, vec![-1, 0, -1]);

// Boundary: node 6 (neutral)
topo.add_node(6, vec![0, 0, 0]);

topo.compute_edges();

// Positive edges: within each cluster
// Negative edges: between clusters
// Zero edges: boundary node is orthogonal to everyone
println!("Positive: {} edges", topo.positive_edges().len());
println!("Negative: {} edges", topo.negative_edges().len());

// Evolve: boundary node stays neutral (no neighbors to pull it)
// Cluster nodes converge toward identical states
let history = topo.evolve(10);
println!("Edge count evolution: {:?}", history);
```

## Design Decisions

**All-to-all affinity** — `compute_edges()` checks every pair of nodes (O(n²)). For small-to-medium networks (hundreds of nodes), this is fine. For large networks, you'd want approximate nearest-neighbor search or locality-sensitive hashing on the ternary state vectors.

**Positive-only neighbor influence** — During evolution, only positive-affinity neighbors pull on a node. Negative-affinity neighbors don't push. This is intentional: the negative edges mark boundaries (opposition), and nodes on the boundary of a cluster shouldn't be pushed away — they should be pulled toward their own cluster. If you want repulsion, add it as a separate force.

**Fixed dimension** — All nodes share the same state vector dimension, set at topology creation. `add_node` enforces this by panicking if a state vector's length differs from the topology dimension, which is what guarantees the affinity function is well-defined for every pair. Variable-dimension nodes would require padding or projection.

**HashMap nodes** — Nodes are stored in a `HashMap<u32, Node>` for O(1) lookup by ID. The IDs are user-assigned (not sequential), which supports dynamic node addition without reindexing.

## Ecosystem Connections

- **`ternary-graph`** — Graph algorithms for ternary-weighted edges (pathfinding, traversal)
- **`ternary-network`** — Network science measures (clustering coefficient, degree distribution)
- **`ternary-geometry`** — Geometric operations on ternary coordinates
- **`ternary-mesh`** — Mesh network topology and routing
- **`ternary-som`** — Self-organizing maps with ternary weights (similar self-organizing principle)
- **`ternary-lattice`** — Lattice structures in ternary spaces
- **`ternary-homology`** — Topological invariants (Betti numbers, Euler characteristic)
- **`ternary-pagerank`** — Centrality measures on ternary graphs

## Open Questions

- **Convergence proof**: Does the evolution always converge? The finite state space ({-1,0,+1}^d per node, finite nodes) guarantees no divergence, but oscillation is possible. Empirically it converges, but a formal proof would be nice.
- **Optimal threshold**: The 0.3 ternarization threshold is ad hoc. What's the principled choice? It might depend on the dimension and network size.
- **Dynamic nodes**: Currently nodes can be added but not removed. A `remove_node()` method would enable dynamic network reconfiguration.
- **Weighted evolution**: Should negative-affinity neighbors exert repulsive force during evolution? The current design ignores them, but opposition-based evolution could create more interesting dynamics.
- **Approximate affinity**: For large networks, computing all n² affinities is expensive. Can locality-sensitive hashing (LSH) on ternary vectors accelerate this?

## Stats

| Metric | Value |
|--------|-------|
| Lines of Rust | ~155 |
| Tests | 11 |
| Dependencies | 0 (uses `std::collections::HashMap`) |

## License

MIT
