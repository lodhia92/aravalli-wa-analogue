"""Strengthen the pathfinder validation.
 (1) Per-element and per-mineral Spearman correlation across pairs (India vs Australia,
     within-survey-standardised), with permutation p-values + Benjamini-Hochberg FDR.
 (2) Composite correlation as the pair set grows: robust-20, all-MNN, all-candidates, and at
     relaxed Australian drainage thresholds (25/15/10%), to show whether the trend strengthens with n.
 Matching is unchanged from analogue_pairing.py (six within-survey-standardised provenance ratios).

Author: Bhavik Harish Lodhia, Curtin University
"""
import os, csv, numpy as np, pandas as pd

import paths

from aravalli_wa.composition import ratios

from aravalli_wa.stats import zscore_elem
from aravalli_wa.constants import P_MASS_FRACTION_OF_P2O5, TI_MASS_FRACTION_OF_TIO2, WT_PCT_TO_MG_KG

def india(name):
    d=pd.read_csv(f"{DR}/{name}_contained.csv"); d=d[d.pct_in_domain>=IN_THR]
    d["k"]=d.lat.round(4).astype(str)+"_"+d.lon.round(4).astype(str)
    m=ar[ar.k.isin(set(d.k))].copy(); m["sid"]=name+"_"+m.k; return m

def findcol(el):
    for meth in ("ICP-MS","XRF"):
        for i,c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"): return i

def aus(path,label,thr):
    d=pd.read_csv(path); d=d[d.pct_in_domain>=thr]
    ids=set(pd.to_numeric(d["id"],errors="coerce").dropna())
    m=ng[ng.index.isin(ids)].copy(); m["sid"]=[f"{label}_{i}" for i in m.index]; return m.reset_index()

def zlog_match(df,mu,sd):
    r=ratios(df).replace([np.inf,-np.inf],np.nan); x=np.log(r.where(r>0)).replace([np.inf,-np.inf],np.nan)
    return (x-mu)/sd

def build_pairs(aus_thr):
    wapp=aus(f"{DR}/wa_palaeoprot_contained.csv","WA_PP",aus_thr)
    yiln=aus(f"{DR}/yilgarn_contained.csv","Y+N",aus_thr)
    ngcm_all=pd.concat([sand,mang]); ngsa_all=pd.concat([wapp,yiln])
    rin=np.log(ratios(ngcm_all).replace([np.inf,-np.inf],np.nan).where(lambda z:z>0))
    rau=np.log(ratios(ngsa_all).replace([np.inf,-np.inf],np.nan).where(lambda z:z>0))
    mi,si=rin.mean(),rin.std(ddof=0); ma,sa=rau.mean(),rau.std(ddof=0)
    out=[]
    for idf,adf,dom in [(sand,wapp,"Palaeoproterozoic"),(mang,yiln,"Archaean")]:
        zi=zlog_match(idf,mi,si)[MATCH].dropna(); za=zlog_match(adf,ma,sa)[MATCH].dropna()
        Iv=zi.values; Av=za.values
        for a in range(len(Av)):
            di=np.sqrt(((Iv-Av[a])**2).sum(1)); bi=int(di.argmin())
            # MNN
            da=np.sqrt(((Av-Iv[bi])**2).sum(1)); mnn=(int(da.argmin())==a)
            out.append(dict(domain=dom,india_sid=idf.loc[zi.index[bi],"sid"],
                            aus_sid=adf.loc[za.index[a],"sid"],dist=float(di[bi]),mnn=mnn))
    P=pd.DataFrame(out)
    return P,ngcm_all,ngsa_all

def spearman(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b); a,b=a[ok],b[ok]
    if len(a)<4: return np.nan,np.nan,len(a)
    ra=pd.Series(a).rank().values; rb=pd.Series(b).rank().values
    rho=float(np.corrcoef(ra,rb)[0,1])
    B=rb[rng.random((NPERM,len(rb))).argsort(axis=1)]
    Bc=B-B.mean(1,keepdims=True); ac=ra-ra.mean()
    den=np.sqrt((Bc**2).sum(1)*(ac**2).sum())
    with np.errstate(invalid="ignore",divide="ignore"):
        null=np.abs(np.where(den>0,(Bc@ac)/den,np.nan))
    p=(np.sum(null>=abs(rho))+1)/(NPERM+1); return rho,p,len(a)

def elem_vecs(P,el,grp=None):
    xs,ys=[],[]
    for _,r in P.iterrows():
        if r.india_sid not in ic.index or r.aus_sid not in acn.index: continue
        di=ic.loc[r.india_sid]; da=acn.loc[r.aus_sid]
        if grp:
            xi=np.nanmean([(np.log(di[e])-mu_in[e])/sd_in[e] for e in grp])
            yi=np.nanmean([(np.log(da[e])-mu_au[e])/sd_au[e] for e in grp])
        else:
            xi=(np.log(di[el])-mu_in[el])/sd_in[el]; yi=(np.log(da[el])-mu_au[el])/sd_au[el]
        xs.append(xi); ys.append(yi)
    return xs,ys

def composite(P):
    xs,ys=[],[]
    for _,r in P.iterrows():
        if r.india_sid not in ic.index or r.aus_sid not in acn.index: continue
        di=ic.loc[r.india_sid]; da=acn.loc[r.aus_sid]
        xs.append(np.nanmean([(np.log(di[e])-mu_in[e])/sd_in[e] for e in PATH]))
        ys.append(np.nanmean([(np.log(da[e])-mu_au[e])/sd_au[e] for e in PATH]))
    return spearman(xs,ys)

def main():
    global CI, DR, FIG, H, HERE, IN_THR, MATCH, MIN, NG, NG_ELEMS, NPERM, P, P25, PATH, PROJ, Pt, RES, SEED, _, acn, ar, c, cidx, crow, d, el, f, grp, ic, k, mang, mu_au, mu_in, n, name, ng, ngcm_all, ngsa_all, ordered, p, rho, rng, robust, rows, sand, sd_au, sd_in, sets, tab, thr, use, v, x, y
    HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(HERE)
    RES=f"{HERE}/results"; DR=f"{RES}/drainage"; FIG=f"{HERE}/figures"
    IN_THR=50.0
    CI=dict(La=.237,Yb=.170,Sm=.148,Eu=.0580,Gd=.199)
    MATCH=["Th/Sc","La/Sc","Th/Co","EuEu","La/Yb_n","Nb/Y"]
    MIN={"zircon (Zr,Hf)":["Zr","Hf"],"monazite (Ce,Nd,Pr)":["Ce","Nd","Pr"],
         "Ti-oxide (Ti)":["Ti"],"xenotime (Dy)":["Dy"],"apatite (P)":["P"]}
    PATH=["Zr","Hf","Ti","Ce","Nd","Pr","Dy","P"]
    NPERM = 100000
    SEED = 20260827
    rng=np.random.default_rng(SEED)
    ar=pd.read_csv(paths.NGCM_TABLE)
    for c in ["Th","Sc","Co","La","Eu","Sm","Gd","Yb","Nb","Y","TiO2","P2O5","Zr","Hf","Ce","Nd","Pr","Dy","LAT","LON"]:
        ar[c]=pd.to_numeric(ar[c],errors="coerce")
    ar["Ti"]=ar["TiO2"]* WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2; ar["P"]=ar["P2O5"]* WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
    ar["k"]=ar.LAT.round(4).astype(str)+"_"+ar.LON.round(4).astype(str)
    sand=india("sandmata"); mang=india("mangalwar")
    NG=paths.NGSA
    with open(NG,encoding="latin-1") as f: H=list(csv.reader(f))[11]
    NG_ELEMS=["Th","Sc","Nb","Y","La","Yb","Co","Eu","Sm","Gd"]+PATH
    cidx={el:findcol(el) for el in NG_ELEMS}; cidx={k:v for k,v in cidx.items() if v is not None}
    ordered=sorted(cidx.items(),key=lambda kv:kv[1]); use=[0]+[v for _,v in ordered]
    ng=pd.read_csv(NG,header=None,skiprows=12,usecols=use,encoding="latin-1",low_memory=False)
    ng.columns=["SITEID"]+[k for k,_ in ordered]
    for c in ng.columns: ng[c]=pd.to_numeric(ng[c],errors="coerce")
    ng=ng.groupby("SITEID").median(numeric_only=True)
    P25,ngcm_all,ngsa_all=build_pairs(25.0)
    mu_in,sd_in=zscore_elem(ngcm_all,PATH); mu_au,sd_au=zscore_elem(ngsa_all,PATH)
    ic=ngcm_all.set_index("sid"); acn=ngsa_all.set_index("sid")
    robust=pd.concat([P25[P25.domain==d].sort_values(["mnn","dist"],ascending=[False,True]).head(10) for d in ["Palaeoproterozoic","Archaean"]])
    print("==== (1) PER-ELEMENT validation on robust-20 ====")
    rows=[]
    for el in PATH:
        x,y=elem_vecs(robust,el); rho,p,n=spearman(x,y); rows.append(("elem",el,rho,p,n))
    for name,grp in MIN.items():
        x,y=elem_vecs(robust,None,grp); rho,p,n=spearman(x,y); rows.append(("mineral",name,rho,p,n))
    tab=pd.DataFrame(rows,columns=["kind","feature","rho","p","n"])
    tab.to_csv(f"{RES}/pair_per_element_validation.csv",index=False)
    print(tab.round(3).to_string(index=False))
    print("\n==== (2) COMPOSITE validation as n grows ====")
    sets=[("robust-20",robust),("all-MNN-25%",P25[P25.mnn]),("all-candidates-25%",P25)]
    for thr in [15.0,10.0]:
        Pt,_,_=build_pairs(thr); sets.append((f"all-candidates-{int(thr)}%",Pt))
    crow=[]
    for name,P in sets:
        rho,p,n=composite(P); crow.append((name,n,round(rho,3),round(p,3)))
        print(f"  {name:22} n={n:3d}  rho={rho:+.3f}  p={p:.3f}")
    pd.DataFrame(crow,columns=["pair_set","n","rho","p"]).to_csv(f"{RES}/pair_composite_scaling.csv",index=False)
    print("\nwrote pair_per_element_validation.csv, pair_composite_scaling.csv")

if __name__ == "__main__":
    main()
