"""Renders a Foehncraft VTKHDF mesh file to a PNG picture for the gallery,
with VTK (the library under ParaView), off screen.

    pip install vtk
    python3 -m tools.mesh.render_mesh FILE.vtkhdf OUT.png [--field NAME]
        [--cut] [--top] [--title TEXT] [--size 900x700] [--range LO HI]

--cut   remove the octant x > xc, y < yc, z > zc to show the inside
--top   look along -z (2D meshes, one cell thick)
--half  remove the cells above the mid plane z = zc (use with --top)
--crop X0 X1 Y0 Y1   keep only the cells inside this x-y window
--no-edges  do not draw the cell edges (fine meshes, flow fields)
--diverging blue-white-red colour map, white at the middle of the range
            (signed fields such as vorticity; use a symmetric --range)
--zoom F    zoom factor of the camera (default 1.25 with --top, else 1.1;
            long thin domains such as a channel need more)
"""
import argparse
import vtk

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("out")
ap.add_argument("--field", default=None)
ap.add_argument("--cut", action="store_true")
ap.add_argument("--top", action="store_true")
ap.add_argument("--half", action="store_true")
ap.add_argument("--crop", nargs=4, type=float, default=None)
ap.add_argument("--title", default="")
ap.add_argument("--size", default="900x700")
ap.add_argument("--range", nargs=2, type=float, default=None)
ap.add_argument("--categorical", action="store_true")
ap.add_argument("--no-edges", action="store_true")
ap.add_argument("--zoom", type=float, default=None)
ap.add_argument("--diverging", action="store_true")
a = ap.parse_args()

r = vtk.vtkHDFReader()
r.SetFileName(a.file)
r.Update()
out = r.GetOutputDataObject(0)
names = [out.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME()) for i in range(out.GetNumberOfPartitionedDataSets())]
grid = out.GetPartitionedDataSet(names.index("internal")).GetPartition(0)
b = grid.GetBounds()
src = grid
if a.cut:
    box = vtk.vtkBox()
    xc, yc, zc = (b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2
    box.SetBounds(xc, b[1] + 1, b[2] - 1, yc, zc, b[5] + 1)
    ex = vtk.vtkExtractGeometry()
    ex.SetInputData(grid)
    ex.SetImplicitFunction(box)
    ex.ExtractInsideOff()
    ex.ExtractBoundaryCellsOff()
    ex.Update()
    src = ex.GetOutput()
def remove_box(data, bounds, inside=False):
    box = vtk.vtkBox()
    box.SetBounds(*bounds)
    ex = vtk.vtkExtractGeometry()
    ex.SetInputData(data)
    ex.SetImplicitFunction(box)
    if inside:
        ex.ExtractInsideOn()
    else:
        ex.ExtractInsideOff()
    ex.ExtractBoundaryCellsOff()
    ex.Update()
    return ex.GetOutput()


if a.half:
    zc = (b[4] + b[5]) / 2
    src = remove_box(src, (b[0] - 1, b[1] + 1, b[2] - 1, b[3] + 1, zc, b[5] + 1))
if a.crop:
    x0, x1, y0, y1 = a.crop
    src = remove_box(src, (x0, x1, y0, y1, b[4] - 1, b[5] + 1), inside=True)
surf = vtk.vtkDataSetSurfaceFilter()
surf.SetInputData(src)
surf.SetNonlinearSubdivisionLevel(0)
surf.Update()

ren = vtk.vtkRenderer()
ren.SetBackground(1, 1, 1)
m = vtk.vtkPolyDataMapper()
m.SetInputConnection(surf.GetOutputPort())
if a.field:
    cd = surf.GetOutput().GetCellData()
    arr = cd.GetArray(a.field)
    if arr is None:
        names = sorted(cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays()))
        raise SystemExit(f"no cell field {a.field!r} in {a.file}; available: {', '.join(names)}")
    m.SetScalarModeToUseCellFieldData()
    m.SelectColorArray(a.field)
    lo, hi = a.range if a.range else arr.GetRange()
    m.SetScalarRange(lo, hi)
    lut = vtk.vtkLookupTable()
    if a.categorical:
        colours = [(0.0, 0.45, 0.70), (0.90, 0.62, 0.0), (0.0, 0.62, 0.45), (0.80, 0.47, 0.65),
                   (0.34, 0.71, 0.91), (0.84, 0.37, 0.0), (0.94, 0.89, 0.26), (0.5, 0.5, 0.5)]
        n = int(hi - lo) + 1
        lut.SetNumberOfTableValues(n)
        for i in range(n):
            lut.SetTableValue(i, *colours[i % len(colours)], 1.0)
        lut.SetRange(lo, hi)
    elif a.diverging:
        ctf = vtk.vtkColorTransferFunction()
        ctf.SetColorSpaceToDiverging()
        ctf.AddRGBPoint(0.0, 0.23, 0.30, 0.75)
        ctf.AddRGBPoint(0.5, 0.87, 0.87, 0.87)
        ctf.AddRGBPoint(1.0, 0.71, 0.02, 0.15)
        lut.SetNumberOfTableValues(256)
        for i in range(256):
            lut.SetTableValue(i, *ctf.GetColor(i / 255.0), 1.0)
        lut.SetRange(lo, hi)
    else:
        lut.SetHueRange(0.667, 0.0)
        lut.SetRange(lo, hi)
    if not a.diverging:
        lut.Build()
    m.SetLookupTable(lut)
    m.ScalarVisibilityOn()
    if not a.categorical:
        bar = vtk.vtkScalarBarActor()
        bar.SetLookupTable(lut)
        bar.SetTitle(a.field.replace("_", " "))
        bar.SetNumberOfLabels(5)
        bar.SetWidth(0.08)
        bar.SetHeight(0.5)
        bar.SetPosition(0.9, 0.25)
        for t in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
            t.SetColor(0.15, 0.15, 0.15)
            t.ShadowOff()
            t.ItalicOff()
            t.BoldOff()
        bar.GetTitleTextProperty().SetFontSize(14)
        ren.AddViewProp(bar)
else:
    m.ScalarVisibilityOff()
act = vtk.vtkActor()
act.SetMapper(m)
if not a.field:
    act.GetProperty().SetColor(0.80, 0.87, 0.97)
if not a.no_edges:
    act.GetProperty().EdgeVisibilityOn()
act.GetProperty().SetEdgeColor(0.12, 0.12, 0.16)
act.GetProperty().SetLineWidth(1.0)
ren.AddActor(act)
if a.title:
    t = vtk.vtkTextActor()
    t.SetInput(a.title)
    t.GetTextProperty().SetColor(0.1, 0.1, 0.1)
    t.GetTextProperty().SetFontSize(20)
    t.SetPosition(15, 15)
    ren.AddViewProp(t)
cam = ren.GetActiveCamera()
if a.top:
    cam.SetPosition(0, 0, 1)
    cam.SetViewUp(0, 1, 0)
else:
    cam.SetPosition(1.6, -1.3, 1.2)
    cam.SetViewUp(0, 0, 1)
cam.SetFocalPoint(0, 0, 0)
ren.ResetCamera()
cam.Zoom(a.zoom if a.zoom else (1.25 if a.top else 1.1))
w, h = (int(x) for x in a.size.split("x"))
win = vtk.vtkRenderWindow()
win.SetOffScreenRendering(1)
win.AddRenderer(ren)
win.SetSize(w, h)
win.Render()
img = vtk.vtkWindowToImageFilter()
img.SetInput(win)
img.Update()
png = vtk.vtkPNGWriter()
png.SetFileName(a.out)
png.SetInputConnection(img.GetOutputPort())
png.Write()
print("wrote", a.out)
