"""Is FMDA's overnight moisture recovery systematically stronger than WRF-SFIRE's?

Finley showed FMDA wetter at every hour, peaking 38% above WRF's own scheme overnight
and nearly converging when dry -- an amplitude error, not a phase one.  If that holds
across the batch it is a bigger lever on 8-hour skill than any pSAF value, because the
extinction threshold turns a small moisture offset into hours of release timing.

Compares WRF-SFIRE's FMC_GC 1-hour class against the FMDA series the forecast actually
used, at the same valid times.  Sampling differs slightly -- WRF is a domain mean over
~30 km, FMDA a radius_km=10 mean about the ignition -- so treat small offsets with care;
the question here is the sign and whether it grows with wetness.
"""
from __future__ import print_function
import sys, glob, os.path as osp, re, warnings
warnings.filterwarnings('ignore'); sys.path[:0]=['src','src/ingest']
import numpy as np, netCDF4 as nc, pandas as pd
import forefire as ff

rows=[]
for line in open(sys.argv[1]):
    parts=line.split('\t')
    if len(parts)<2: continue
    job, iss = parts[0].strip(), parts[1].strip()
    w=osp.dirname(job)
    try:
        _,ign,latlon,_=ff.read_input(w)
        tt=ff.make_timing_table(w,ign)
        ser=ff.md_series(tt,latlon,ff.load_ff_cfg(),issue_utc=iss)
    except Exception:
        continue
    if ser is None: continue
    fm={}
    for i,(_,r) in enumerate(tt.iterrows()):
        if r['ign_seconds']<0 or ser[i] is None: continue
        fm[pd.Timestamp(str(r['UTC_str']).replace('Z',''))]=ser[i][0]
    outs=sorted(glob.glob(osp.join(w,'wrf','wrfout_d01_*'))) or \
         sorted(glob.glob(osp.join(w,'wrf','saveout*')))
    pairs=[]
    for p in outs[::2]:
        m=re.search(r'_(\d{4}-\d\d-\d\d)_(\d\d):(\d\d)',osp.basename(p))
        if not m: continue
        t=pd.Timestamp('%sT%s:%s:00'%(m.group(1),m.group(2),m.group(3)))
        if t not in fm: continue
        try:
            d=nc.Dataset(p)
            if 'FMC_GC' not in d.variables: d.close(); continue
            a=np.array(d['FMC_GC'][:]).squeeze(); d.close()
        except Exception: continue
        wv = float(a[0].mean()) if a.shape[0]<=6 else float(a[:,:,0].mean())
        pairs.append((t,wv,fm[t]))
    if len(pairs)<6: continue
    W=np.array([p[1] for p in pairs]); F=np.array([p[2] for p in pairs])
    dry=W<=np.percentile(W,25); wet=W>=np.percentile(W,75)
    rows.append((osp.basename(w)[4:].split('_20')[0], len(pairs),
                 (F-W).mean(), (F-W)[wet].mean(), (F-W)[dry].mean(),
                 W.max()-W.min(), F.max()-F.min()))

print('%-26s %4s %9s %10s %10s %10s %10s'
      %('fire','n','mean diff','when wet','when dry','WRF range','FMDA range'))
for r in sorted(rows,key=lambda x:-x[3]):
    print('%-26s %4d %+9.4f %+10.4f %+10.4f %10.4f %10.4f'%r)
A=np.array([[r[2],r[3],r[4],r[5],r[6]] for r in rows])
print('\n%-26s %4d %+9.4f %+10.4f %+10.4f %10.4f %10.4f'
      %('MEAN',len(rows),A[:,0].mean(),A[:,1].mean(),A[:,2].mean(),A[:,3].mean(),A[:,4].mean()))
print('\nFMDA wetter overall : %d of %d fires'%((A[:,0]>0).sum(),len(A)))
print('wet-end gap exceeds dry-end gap: %d of %d'%((A[:,1]>A[:,2]).sum(),len(A)))
print('FMDA diurnal swing larger       : %d of %d'%((A[:,4]>A[:,3]).sum(),len(A)))
