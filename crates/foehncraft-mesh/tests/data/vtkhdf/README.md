# VTKHDF fixtures

Written by KoflCraft (PR #38 of chiefenne/KoflCraft, the blocked VTKHDF
writer), an independent writer of the blocked layout: VTKHDF 2.4,
`PartitionedDataSetCollection`, block `internal` with one partition per
block of owned cells plus a ghost layer, one `PolyData` block per patch.

- `box4.vtkhdf`, `box4.vtu`: KoflCraft's box-4 fixture exported both ways
  (4 cells, 6 patches, the `koflcraft_*` digests as root attributes).
- `channel.vtkhdf`, `channel.vtu`, `channel-mesh-evidence.json`: a straight
  channel of 776 cells in one block, patches `inlet`, `outlet` and `walls`;
  the evidence file lists the SHA-256 of `channel.vtkhdf` (KoflCraft calls
  it `foehncraft.vtkhdf`) and the summary the reader must reproduce.

Read by `tests/vtkhdf_read.rs` and the corpus of `tests/broken_inputs.rs`.
