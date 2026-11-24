"""Draws the images of the stirred-tank entry from a finished run: speed in
the vertical plane y = 0 with streamlines, speed in a horizontal plane below
the impeller, and the pressure on the impeller.

    pip install h5py numpy scipy matplotlib
    python3 examples/applications/stirred_tank/images.py output/case docs/figures/applications

The run directory needs `fields.vtkhdf` (output.vtkhdf = true); the tank is
T = 1 m, 90 rpm (mesh.py).
"""
import sys, h5py, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))
from foehncraft.vtkhdf import internal  # noqa: E402  (fields.vtkhdf in one piece)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3]))
from tools.plot.colormaps import RAINBOW  # noqa: E402
run, out = sys.argv[1], sys.argv[2]
TIP = 9.42477796076938*0.2
f = h5py.File(run+"/fields.vtkhdf","r"); g=internal(f)
pts=g["Points"][:]; con=g["Connectivity"][:].reshape(-1,4)
cen=pts[con].mean(axis=1); cd=g["CellData"]
u=np.stack([cd["velocity_x"][:],cd["velocity_y"][:],cd["velocity_z"][:]],1)
tree=cKDTree(cen)

HEADD=0.25; RT=0.5
def inside(R,Z):
    r=np.abs(R); zb=HEADD*(1-np.sqrt(np.clip(1-(r/RT)**2,0,1)))
    return (r<=RT)&(Z>=zb)&((r>0.026)|(Z<0.285))
def grid(plane, n=520):
    if plane=="v":
        r=np.linspace(-0.5,0.5,n); z=np.linspace(0,1.0,n); R,Z=np.meshgrid(r,z)
        P=np.stack([R,0*R,Z],-1).reshape(-1,3); msk=inside(R,Z).reshape(-1)
    d,i=tree.query(P); ok=msk&(d<0.06)
    return R,Z,i,ok
R,Z,i,ok=grid("v")
U=u[i]; sp=np.linalg.norm(U,axis=1)
mm=lambda v: np.where(ok,v,np.nan).reshape(R.shape)
fig,ax=plt.subplots(figsize=(6.2,6.4),dpi=200)
im=ax.pcolormesh(R,Z,mm(sp)/TIP,cmap=RAINBOW,vmin=0,vmax=0.5,shading="auto",rasterized=True)
ax.streamplot(R[0],Z[:,0],np.nan_to_num(mm(U[:,0])),np.nan_to_num(mm(U[:,2])),color="k",linewidth=0.5,density=2.0,arrowsize=0.6)
ax.set_aspect("equal"); ax.set_xlabel("x  (m)"); ax.set_ylabel("z  (m)")
cb=fig.colorbar(im,ax=ax,shrink=0.8,pad=0.02); cb.set_label("speed / tip speed  (-)")
ax.set_title("Speed and streamlines in the vertical plane y = 0",fontsize=10)
fig.tight_layout(); fig.savefig(out+"/stirred_tank_speed_section.png"); plt.close(fig)
# horizontal section below the impeller
zs=0.22; n=520; x=np.linspace(-0.5,0.5,n); X,Y=np.meshgrid(x,x)
P=np.stack([X,Y,zs+0*X],-1).reshape(-1,3); d,i=tree.query(P); ok=((X**2+Y**2)<=0.25).reshape(-1)&(d<0.06)
U=u[i]; sp=np.linalg.norm(U,axis=1)
fig,ax=plt.subplots(figsize=(6.2,5.6),dpi=200)
im=ax.pcolormesh(X,Y,np.where(ok,sp,np.nan).reshape(X.shape)/TIP,cmap=RAINBOW,vmin=0,vmax=0.5,shading="auto",rasterized=True)
th=np.linspace(0,2*np.pi,300); ax.plot(0.5*np.cos(th),0.5*np.sin(th),"k",lw=0.8)
ax.set_aspect("equal"); ax.set_xlabel("x  (m)"); ax.set_ylabel("y  (m)")
cb=fig.colorbar(im,ax=ax,shrink=0.8,pad=0.02); cb.set_label("speed / tip speed  (-)")
ax.set_title("Speed in the horizontal plane z = 0.22 m, below the impeller",fontsize=10)
fig.tight_layout(); fig.savefig(out+"/stirred_tank_speed_horizontal.png"); plt.close(fig)
# impeller close-up coloured by pressure
def surf(name):
    s=f["VTKHDF/"+name]; P=s["Points"][:]; c=s["Polygons/Connectivity"][:]; o=s["Polygons/Offsets"][:]
    return [P[c[o[k]:o[k+1]]] for k in range(len(o)-1)], s["CellData"]
fig=plt.figure(figsize=(6.5,6),dpi=200); ax=fig.add_subplot(111,projection="3d")
polys,cdi=surf("impeller")
cen_i=np.array([p.mean(0) for p in polys])
# pressure from nearest fluid cell
pr=cd["pressure"][:]; _,ii=tree.query(cen_i); pv=pr[ii]
lo,hi=np.percentile(pv,[3,97])
ax.add_collection3d(Poly3DCollection(polys,facecolors=plt.cm.coolwarm((pv-lo)/(hi-lo)).clip(0,1),edgecolors="none"))
for name in ("shaft","baffles","sparger"):
    pl,_=surf(name); ax.add_collection3d(Poly3DCollection(pl,facecolors="0.72",edgecolors="none",alpha=0.9))
ax.set_xlim(-0.25,0.25); ax.set_ylim(-0.25,0.25); ax.set_zlim(0.1,0.6); ax.set_box_aspect((1,1,1)); ax.view_init(elev=24,azim=-60); ax.set_axis_off()
sm=plt.cm.ScalarMappable(cmap="coolwarm",norm=plt.Normalize(lo,hi)); cb=fig.colorbar(sm,ax=ax,shrink=0.6,pad=0.0); cb.set_label("pressure on the impeller  (Pa)")
ax.set_title("Three-blade hydrofoil impeller, hub, shaft and sparger ring",fontsize=10)
fig.savefig(out+"/stirred_tank_impeller.png"); plt.close(fig)
