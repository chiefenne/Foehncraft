- **Faster CI.** The integration tests of `foehncraft-flow` are two test
  programs, `flow` and `flow2`, that compile side by side (the single
  program compiled for nine minutes on one thread of the compiler front
  end); the test job runs all test programs in one pool with
  cargo-nextest; the MPI job builds its example programs in one build with
  the `ci` profile instead of eight `release` builds.
