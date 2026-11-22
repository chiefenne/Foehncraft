"""Images of the heat-exchanger mesh: a section through all six zones and a
cut-away of the housing with the tube bundle (pyvista, off-screen)."""
import sys, numpy as np, meshio, pyvista as pv
pv.OFF_SCREEN=True
m=meshio.read(sys.argv[1]); outdir=sys.argv[2]
names={v[0]:k for k,v in m.field_data.items() if v[1]==3}
cells=[];reg=[]
for blk,phys in zip(m.cells,m.cell_data["gmsh:physical"]):
    if blk.type=="tetra": cells.append(blk.data); reg.append(phys)
g=pv.UnstructuredGrid({pv.CellType.TETRA:np.vstack(cells)},m.points*1000); g.cell_data["region"]=np.concatenate(reg)
col={"shell_side":"#f4a261","tube_side":"#4ea8de","tubes":"#c0392b","tubesheets":"#e9c46a","baffles":"#d4a017","housing":"#8d99ae"}
lab={"shell_side":"shell side (hot air)","tube_side":"tube side (water)","tubes":"tubes","tubesheets":"tubesheets","baffles":"baffles","housing":"housing"}
def zone(grid,n):
    tag=[k for k,v in names.items() if v==n][0]
    return grid.extract_cells(np.where(grid.cell_data["region"]==tag)[0])
cam=[(470,-190,300),(0,150,25),(0,0,1)]
def plotter():
    p=pv.Plotter(window_size=(2400,1200)); p.set_background("#1b1f2a"); p.enable_anti_aliasing("ssaa"); return p
# 1. section: every cell with its centre at x < 0 (whole cells, the mesh shows)
half=g.extract_cells(np.where(g.cell_centers().points[:,0]<0)[0])
p=plotter()
for n in ["housing","shell_side","tube_side","tubes","tubesheets","baffles"]:
    p.add_mesh(zone(half,n).extract_surface(algorithm="dataset_surface"),color=col[n],show_edges=True,edge_color="#20242c",line_width=0.4,label=lab[n])
p.add_legend(bcolor="#1b1f2a",border=False,size=(0.2,0.24),loc="lower right",face="rectangle")
p.camera_position=cam; p.camera.zoom(1.15); p.screenshot(f"{outdir}/hx_mesh_section.png")
# 2. cut-away: half of the housing, the whole tube bundle with its surface mesh
p=plotter()
p.add_mesh(zone(half,"housing").extract_surface(algorithm="dataset_surface"),color=col["housing"],show_edges=True,edge_color="#3a404c",line_width=0.4)
for n in ["tubes","tubesheets","baffles"]:
    p.add_mesh(zone(g,n).extract_surface(algorithm="dataset_surface"),color=col[n],show_edges=True,edge_color="#2a1a14",line_width=0.4,specular=0.3)
p.camera_position=cam; p.camera.zoom(1.15); p.screenshot(f"{outdir}/hx_mesh_cutaway.png")
