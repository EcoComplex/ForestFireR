# NEWS

## ForestFireR 0.4.01

### Bug fixes

* Fixed `InitialConditionNonHomogeneous()` (used by `generate_landscape_layers()`
  for multi-species initial landscapes, e.g. Fig. 8 and Fig. 9 of the
  supplementary material). Pair classification during the clustering
  growth process relied on numeric-mean thresholds
  (`threshold1 = (state1+state1)/2`, etc.) that only classify pairs
  correctly when `state1 < state2`. For native/invader layers laid over an
  `empty_postfire` background (`background = 5 > pattern states 1, 2`),
  this assumption did not hold, corrupting cluster counts for every
  clustering parameter `p` and silently producing landscapes with the
  wrong spatial structure (e.g. very different `p` values could produce
  nearly identical patterns at `t = 0`).
* Replaced the threshold-based classification with direct state-equality
  comparisons (`S[i][j] == state1`, etc.) in both the pair-counting loop
  and the two move-evaluation blocks (`p1`, `p2`) of
  `InitialConditionNonHomogeneous()`. Verified by cluster-count regression
  tests against the single-layer reference implementation, and by a full
  paper-scale re-run of Fig. 9 showing the expected monotonic ordering by
  `p` and the expected cluster shapes.
* Only initial-landscape generation for **two or more layers** was
  affected (Fig. 8 and Fig. 9's reduced-dynamics scenarios); single-layer
  landscapes (Fig. 10, Fig. 11) were not affected.

### Investigated, not a bug

* At `eta_inv = 0` (fully symmetric parameters, no fire, no replacement),
  `rho_nat` for `p = 1` sits well below the ~0.5 expected by symmetry at
  `T = 300` (Fig. 8/9's simulation horizon), while lower-`p` curves are
  already close to 0.5 -- visible in the original published Fig. 8 too.
  Investigated and confirmed this is NOT an asymmetry between native and
  invader: at `T = 300`, `p = 1` gives `n1 ~= n2 ~= 0.30` with ~40% of the
  grid still unconverted empty-post-fire background (`n5 ~= 0.40`);
  letting the same runs continue to `T ~= 1000` brings both `n1` and `n2`
  to ~0.50 with `n5 = 0`. Root cause: colonization of empty-post-fire
  cells is driven by local perimeter (neighbour count), and a fully
  clustered (`p = 1`) landscape has far less perimeter per occupied cell
  than a fragmented one, so it consumes the remaining background much
  more slowly. `T = 300` is simply not long enough for the `p = 1` curve
  to reach the steady state that lower-`p` curves already reached --
  native and invader are affected equally (`n1` and `n2` track each other
  throughout). Not a code defect; a modeling/figure-design consideration
  (whether Fig. 8/9's fixed `T` should be treated as a snapshot of
  clustering-dependent transient dynamics, or whether the `p = 1` curve
  needs a longer horizon to compare like-for-like against the others).

## ForestFireR 0.4.0

* Initial versioned release used for the manuscript figures (mean-field
  and spatial Gillespie simulation engine, R wrapper functions using the
  paper's own notation, figure-reproduction scripts).
