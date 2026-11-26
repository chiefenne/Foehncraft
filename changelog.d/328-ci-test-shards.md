- **Faster CI.** The tests run in three shards per platform (the flow
  tests in two halves and every other crate) with cargo-nextest, each
  compiling only what it runs, and the MPI job builds its example
  programs in one build with the `ci` profile instead of eight `release`
  builds.
