#!/usr/bin/env python3
"""Recreate the finite-spin threshold plot included in the paper."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify import c3,c4

r=np.linspace(1,10,400)
fig,ax=plt.subplots(figsize=(6.5,3.6))
ax.plot(r,np.ones_like(r),label=r'$s=1$')
ax.plot(r,r/(r+1),label=r'$s=2$')
ax.plot(r,np.array([c3(float(t),1.0) for t in r]),label=r'$s=3$')
ax.plot(r,np.array([c4(float(t),1.0) for t in r]),label=r'$s=4$')
ax.set_xlabel(r'Anisotropy ratio $r=a/b$')
ax.set_ylabel(r'Normalized threshold $c_s/b$')
ax.set_xlim(1,10)
ax.set_ylim(0,1.07)
ax.grid(True,alpha=0.25)
ax.legend(frameon=False,ncol=4,fontsize=9,loc='lower center',bbox_to_anchor=(0.5,1.02))
fig.tight_layout()
fig.savefig('thresholds.pdf',bbox_inches='tight')
fig.savefig('thresholds.png',dpi=180,bbox_inches='tight')
print('Generated thresholds.pdf and thresholds.png')
