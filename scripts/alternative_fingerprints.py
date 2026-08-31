"""
Are broader multi-element rare-earth fingerprints better proxies for their host minerals than
the published sets? Specifically, is a multi-element heavy rare-earth fingerprint a better proxy
for xenotime than Dy alone, and is Ce+Nd+Pr the right light rare-earth set for monazite? Tested against observed HMMA grain proportions
at the drainage-selected Australian sites, same join as mineralogical_validation.py.
Output: results/alternative_fingerprints.csv

Author: Bhavik Harish Lodhia, Curtin University
"""
import csv, os

import numpy as np, pandas as pd

import paths

from aravalli_wa.stats import perm_p

def spearman(x,y):
    xr=pd.Series(x).rank().values; yr=pd.Series(y).rank().values
    xr=xr-xr.mean(); yr=yr-yr.mean()
    d=np.sqrt((xr**2).sum()*(yr**2).sum())
    return float((xr*yr).sum()/d) if d>0 else np.nan

def main():
    global AUS_THR, CAND, DR, ELS, H, HERE, METH, NG, NPERM, PROJ, RES, _, c, cidx, d, e, el, els, f, hm, i, ids, j, k, label, m, missing, ng, o, ordered, pairs, pcol, pid, pool, rho, rows, s, scale, stem, v, x, z
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(os.path.dirname(HERE))
    RES  = os.path.join(HERE, "results"); DR = os.path.join(RES, "drainage")
    NG   = paths.NGSA
    AUS_THR = 25.0; NPERM = 100000
    CAND = {
     "Xenotime: Dy (as published)": (["Dy"], "Xenotime-Y"),
     "Xenotime: Y":                 (["Y"], "Xenotime-Y"),
     "Xenotime: Y+Dy":              (["Y","Dy"], "Xenotime-Y"),
     "Xenotime: Dy+Ho+Er+Yb+Lu":    (["Dy","Ho","Er","Yb","Lu"], "Xenotime-Y"),
     "Xenotime: Y+Dy+Ho+Er+Yb+Lu":  (["Y","Dy","Ho","Er","Yb","Lu"], "Xenotime-Y"),
     "Monazite: Ce+Nd+Pr (as published)": (["Ce","Nd","Pr"], "Monazite"),
     "Monazite: Ce":                (["Ce"], "Monazite"),
     "Monazite: La+Ce+Pr+Nd+Sm":    (["La","Ce","Pr","Nd","Sm"], "Monazite"),
     "Monazite: Ce+Nd+Pr+Th":       (["Ce","Nd","Pr","Th"], "Monazite"),
     "Monazite: Ce+Nd+Pr+Th+P":     (["Ce","Nd","Pr","Th","P"], "Monazite"),
     "Allanite: La+Ce+Ca":          (["La","Ce","Ca"], "Allanite"),
    }
    ELS = sorted({e for v in CAND.values() for e in v[0]})
    METH = {"Ti":"XRF","P":"XRF","Ca":"XRF"}
    with open(NG,encoding="latin-1") as f: H=list(csv.reader(f))[11]
    cidx={}
    for el in ELS:
        m=METH.get(el,"ICP-MS")
        for i,c in enumerate(H):
            if c.strip().startswith("%s %s"%(el,m)): cidx[el]=i; break
    missing=[e for e in ELS if e not in cidx]
    print("elements not found in NGSA:",missing if missing else "none")
    ordered=sorted(cidx.items(),key=lambda kv:kv[1])
    ng=pd.read_csv(NG,header=None,skiprows=12,usecols=[0]+[v for _,v in ordered],
                   encoding="latin-1",low_memory=False)
    ng.columns=["SITEID"]+[k for k,_ in ordered]
    for c in ng.columns: ng[c]=pd.to_numeric(ng[c],errors="coerce")
    ng=ng.groupby("SITEID").median(numeric_only=True).reset_index()
    ids=set()
    for f in ["wa_palaeoprot_contained.csv","yilgarn_contained.csv"]:
        d=pd.read_csv(os.path.join(DR,f))
        ids|=set(pd.to_numeric(d[d.pct_in_domain>=AUS_THR]["id"],errors="coerce").dropna())
    pool=ng[ng.SITEID.isin(ids)].copy()
    z=pd.DataFrame({"SITEID":pool.SITEID.values})
    for el in cidx:
        v=pd.to_numeric(pool[el],errors="coerce"); x=np.log(v.where(v>0))
        z[el]=((x-x.mean())/x.std(ddof=0)).values
    hm=pd.read_csv(os.path.join(RES,"hmma_grain_counts.csv"))
    j=z.merge(hm,on="SITEID",how="inner")
    pairs=pd.read_csv(os.path.join(RES,"analogue_pairs.csv"))
    pairs=pairs[pairs.mnn.astype(str).str.lower()=="true"]
    pid=set(float(str(s).split("_")[-1]) for s in pairs.aus_sid)
    j["is_pair"]=j.SITEID.isin(pid)
    rows=[]
    for label,(els,stem) in CAND.items():
        els=[e for e in els if e in cidx]
        if not els: continue
        pcol="%s (pmo)"%stem
        if pcol not in j.columns: 
            print("no HMMA column for",stem); continue
        j["_s"]=j[els].mean(axis=1)
        for scale,d in [("pool (n=85)",j),("pairs (n=20)",j[j.is_pair])]:
            m=d[["_s",pcol]].apply(pd.to_numeric,errors="coerce").dropna()
            if len(m)<6: continue
            rho=spearman(m["_s"].values,m[pcol].values)
            rows.append(dict(scale=scale,fingerprint=label,elements="+".join(els),
                             n=len(m),rho=round(rho,3),perm_p=round(perm_p(m["_s"].values,m[pcol].values,rho),4)))
    o=pd.DataFrame(rows)
    o.to_csv(os.path.join(RES,"alternative_fingerprints.csv"),index=False)
    print()
    print(o[o.scale=="pool (n=85)"].to_string(index=False))
    print()
    print(o[o.scale=="pairs (n=20)"].to_string(index=False))

if __name__ == "__main__":
    main()
