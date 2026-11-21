# Broken mesh files (D202)

Hand-written inputs that every mesh reader must refuse with an error
message, never a panic, an abort or a hang
(`tests/broken_inputs.rs`, `hand_written_broken_inputs_are_refused`).
The file name starts with the reader (`gmsh`, `su2`, `industry`, `vtk`) and
says what is broken. Most of them reproduce a failure that the generated
corpus of the same test found: a count or index beyond the size of the
file (an allocation that aborted), a header without its newline (a read
past the end), an element block whose type does not fit its dimension.
