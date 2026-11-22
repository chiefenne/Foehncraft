- **Fix: the bounded face velocity without a limiter crashed (D50).**
  `solver.face_velocity = "bounded"` set by hand with no momentum limiter
  read an empty limiter array and ended in a segmentation fault; the
  limiter's per-cell epsilon is now set up whenever the face velocity uses
  it, and an empty array is an error instead of a crash. The automatic
  remediation was not affected.
