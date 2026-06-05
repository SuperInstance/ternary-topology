# ternary-topology

Persistent homology for ternary networks.

Computes Betti numbers for 0-state clusters, builds Vietoris-Rips complexes, detects boundaries between state regions, and identifies topological insulators (0-clusters as barriers separating -1 and +1 regions).

## Features

- **Vietoris-Rips complex**: Build simplicial complexes from point clouds with filtration
- **Betti number computation**: Mod-2 homology via boundary matrix reduction
- **Cluster analysis**: Find connected components by ternary state
- **Boundary detection**: Identify points at state interfaces
- **Topological insulator detection**: Zero-state clusters that separate opposing states
- **Distance-based topology**: All computations driven by spatial proximity

## Usage

```rust
use ternary_topology::{TernaryPoint, Ternary, VietorisRips, ClusterAnalyzer};

let points = vec![
    TernaryPoint::new(0, vec![0.0, 0.0], Ternary::Neg),
    TernaryPoint::new(1, vec![1.0, 0.0], Ternary::Zero),
    TernaryPoint::new(2, vec![2.0, 0.0], Ternary::Pos),
];

// Build Vietoris-Rips complex
let vr = VietorisRips::build(&points, 2.0, 2);

// Compute Betti numbers
let b0 = ClusterAnalyzer::betti_0(&points, 2.0);
let b1 = ClusterAnalyzer::betti_1(&points, 2.0);
println!("Betti-0: {}, Betti-1: {}", b0, b1);

// Detect topological insulators
let insulators = ClusterAnalyzer::detect_insulators(&points, 1.5);
println!("Insulator clusters: {}", insulators.len());
```

## Test Coverage

20 tests covering point distances, simplex operations, VR complex construction, cluster analysis, boundary detection, insulator detection, and Betti number computation.

## Known Limitations

- Betti number computation uses simplified algorithm, not full persistent homology
- Vietoris-Rips complex construction is O(n³) for 2-simplices, not scalable to large datasets
- No barcode or persistence diagram output
- Betti-1 uses Euler characteristic approximation
- Boundary matrix reduction is naive Gaussian elimination, not optimized
- No support for weighted complexes or filtered simplicial complexes beyond distance

## See Also

- **ternary-graph** — Graph algorithms for ternary-weighted edges
- **ternary-network** — Network science for ternary-weighted graphs
- **ternary-geometry** — Geometric operations on ternary coordinates
- **ternary-mesh** — Mesh network topology and routing
- **ternary-som** — Self-organizing maps with ternary weights
- **ternary-lattice** — Lattice structures in ternary spaces

## License

MIT
