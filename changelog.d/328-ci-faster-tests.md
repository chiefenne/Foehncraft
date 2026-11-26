- **Faster CI.** The integration tests of `foehncraft-flow` are two test
  programs, `flow` and `flow2`, that compile side by side (the single
  program compiled for nine minutes on one thread of the compiler front
  end); the tests run with cargo-nextest in one pool, in two halves per
  platform side by side; the three Linux case sets are balanced on their
  run times (the third took 15 minutes, the first 4); the MPI job builds
  its example programs in one build with the `ci` profile instead of eight
  `release` builds.
