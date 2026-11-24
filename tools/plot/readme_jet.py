"""Renders the impinging-jet tile of the README (assets/readme/impinging_jet_temperature.jpg).

    target/release/foehn flow impinging_jet.toml   # with [output] vtkhdf = true
    python3 tools/plot/readme_jet.py output/impinging_jet/fields.vtkhdf jet.png 1.5 split

Left half: log10 of the passive temperature theta; right half: the speed;
white streamlines; the axisymmetric wedge mirrored about the jet axis.
Needs h5py, NumPy, SciPy and Matplotlib.
"""
import sys, h5py, numpy as np, matplotlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
from foehncraft.vtkhdf import internal  # the blocks of fields.vtkhdf in one piece
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import griddata
try:  # tools/plot/colormaps.py, run as a module or as a script
    from tools.plot.colormaps import RAINBOW
except ImportError:
    from colormaps import RAINBOW
f=h5py.File(sys.argv[1]); g=internal(f)
P=g['Points'][:]; C=g['Connectivity'][:]; O=g['Offsets'][:]
cen=np.array([P[C[O[i]:O[i+1]]].mean(0) for i in range(len(O)-1)])
x=cen[:,0]; r=np.hypot(cen[:,1],cen[:,2])
d=g['CellData']; ux=d['velocity_x'][:]; uy=d['velocity_y'][:]; uz=d['velocity_z'][:]
ur=(uy*cen[:,1]+uz*cen[:,2])/np.maximum(r,1e-12); sp=d['speed'][:]; th=d['theta'][:]
R=float(sys.argv[3]); 
# image grid: horizontal = r in [-R,R], vertical = x in [0,2]
nx,ny=900,600
rr=np.linspace(-R,R,nx); xx=np.linspace(0,2,ny); RR,XX=np.meshgrid(rr,xx)
pts=np.c_[x,r]; q=np.c_[XX.ravel(),np.abs(RR.ravel())]
def I(v):
    a=griddata(pts,v,q,method='linear'); b=griddata(pts,v,q,method='nearest'); a[np.isnan(a)]=b[np.isnan(a)]; return a.reshape(ny,nx)
S=I(sp); T=I(th); UX=I(ux); UR=I(ur)*np.sign(RR)
mode=sys.argv[4]
fig=plt.figure(figsize=(6,4),dpi=150); ax=fig.add_axes([0,0,1,1]); bg='#0b0f17'; fig.patch.set_facecolor(bg); ax.set_facecolor(bg)
if mode=='split':
    Tn=np.log10(np.maximum(T,1e-6)); L=RR<0
    ax.imshow(np.where(L,Tn,np.nan),origin='lower',extent=[-R,R,0,2],cmap='inferno',vmin=np.nanpercentile(Tn,8),vmax=Tn.max(),aspect='auto',interpolation='bilinear')
    ax.imshow(np.where(~L,S,np.nan),origin='lower',extent=[-R,R,0,2],cmap=RAINBOW,aspect='auto',interpolation='bilinear')
    ax.plot([0,0],[0.03,2],color=(1,1,1,0.35),lw=0.6,ls='--')
    kw=dict(color='white',fontsize=8,alpha=0.85,va='top')
    ax.text(-R+0.06,1.95,'temperature',**kw); ax.text(R-0.06,1.95,'speed',ha='right',**kw)
elif mode=='theta':
    Tn=np.log10(np.maximum(T,1e-6)); 
    ax.imshow(np.where(np.isnan(Tn),np.nan,Tn),origin='lower',extent=[-R,R,0,2],cmap='inferno',vmin=np.nanpercentile(Tn,40),vmax=np.nanmax(Tn),aspect='auto',interpolation='bilinear')
else:
    ax.imshow(S,origin='lower',extent=[-R,R,0,2],cmap=RAINBOW,aspect='auto',interpolation='bilinear')
ax.streamplot(rr,xx,np.nan_to_num(UR),np.nan_to_num(UX),color=(1,1,1,0.55),density=1.6,linewidth=0.6,arrowsize=0.6)
# plate coloured by wall temperature (heat transfer)
ax.add_patch(plt.Rectangle((-R,-0.0),2*R,0.025,color='#9aa4b5'))
for s in [-1,1]: ax.add_patch(plt.Rectangle((s*0.5 if s>0 else -0.56,1.9),0.06,0.1,color='#9aa4b5'))
ax.set_xlim(-R,R); ax.set_ylim(0,2); ax.axis('off')
fig.savefig(sys.argv[2],facecolor=bg)
