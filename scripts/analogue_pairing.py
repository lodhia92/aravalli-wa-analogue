"""
Robust analogue PAIRS between drainage-defined Indian (NGCM) and Australian (NGSA)
catchments, matched on within-survey-standardised weathering-robust provenance ratios.

Two age domains, matched separately:
  Palaeoproterozoic: Sandmata (India)  <->  WA Palaeoproterozoic (Halls Creek+Lamboo+Capricorn+Gascoyne)
  Archaean:          Mangalwar (India) <->  Youanmi+Narryer (Yilgarn)

Method (documented in the project log, section 5c):
  1. For each matching ratio r, work in natural log: x = ln(r).
  2. Standardise WITHIN each survey s (NGSA, NGCM) to its own mean/SD (Singer & Kouda 2001
     z-score; removes the systematic inter-survey offset that levelling would otherwise correct -
     Grunsky 2010; Grunsky & de Caritat 2020):   z = (x - mean_s) / sd_s.
  3. Distance between India catchment i and Australia catchment j:
       d(i,j) = sqrt( sum_k (z_i,k - z_j,k)^2 )   (Euclidean in standardised log-ratio space).
  4. For each Australian catchment (the smaller pool) take its nearest Indian catchment; flag
     mutual nearest neighbours (MNN) as the most robust pairs; rank by d.
  5. Validate INDEPENDENTLY on critical-mineral pathfinder elements (Zr,Hf zircon; Ti Ti-oxides;
     Ce,Nd,Pr monazite-LREE; Dy xenotime-HREE; P apatite) - elements not used in matching.

Matching ratios (weathering-robust source-rock discriminators; McLennan/Taylor):
  Th/Sc, La/Sc, Th/Co, Eu/Eu*, (La/Yb)n, Nb/Y.
Chondrite (McDonough & Sun 1995, CI): La .237 Yb .170 Sm .148 Eu .0580 Gd .199 ppm.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import os, csv

import numpy as np, pandas as pd

import paths

from aravalli_wa.composition import ratios

def india(name):
    d=pd.read_csv(os.path.join(DR,f"{name}_contained.csv")); d=d[d.pct_in_domain>=IN_THR]
    d["k"]=d.lat.round(4).astype(str)+"_"+d.lon.round(4).astype(str)
    m=ar[ar.k.isin(set(d.k))].copy(); m["sid"]=name+"_"+m.k; m["survey"]="NGCM"
    return m

def col(el,method):
    for i,c in enumerate(H):
        if c.strip().startswith(f"{el} {method}"): return i
    return None

def aus(path,label):
    d=pd.read_csv(path); d=d[d.pct_in_domain>=AUS_THR]
    ids=set(pd.to_numeric(d["id"],errors="coerce").dropna())
    m=ng[ng.SITEID.isin(ids)].copy(); m["sid"]=label+"_"+m.SITEID.astype(str); m["survey"]="NGSA"
    return m

def zlog(df,muSD=None):
    r=ratios(df).replace([np.inf,-np.inf],np.nan)
    x=np.log(r.where(r>0)).replace([np.inf,-np.inf],np.nan)
    if muSD is None: muSD=(x.mean(skipna=True),x.std(ddof=0,skipna=True))
    z=(x-muSD[0])/muSD[1]; return z,muSD

def pair(india_df,aus_df,domain):
    # survey-level mean/SD from the WHOLE drainage pool of each survey (not just this domain)
    zin,_=zlog(india_df, ngcm_muSD); zau,_=zlog(aus_df, ngsa_muSD)
    zin=zin.dropna(); zau=zau.dropna()
    Iv=zin[MATCH].values; Av=zau[MATCH].values
    rows=[]
    for a in range(len(Av)):
        dd=np.sqrt(((Iv-Av[a])**2).sum(1)); bi=int(dd.argmin())
        rows.append(dict(domain=domain,
            aus_sid=aus_df.loc[zau.index[a],"sid"], india_sid=india_df.loc[zin.index[bi],"sid"],
            dist=round(float(dd[bi]),3),
            aus_idx=zau.index[a], india_idx=zin.index[bi]))
    p=pd.DataFrame(rows)
    # mutual nearest neighbour flag: is the India member's nearest Australian the paired one?
    for n,row in p.iterrows():
        iv=zin.loc[row.india_idx,MATCH].values.astype(float)
        dd=np.sqrt(((Av-iv)**2).sum(1)); p.loc[n,"mnn"]=(zau.index[int(dd.argmin())]==row.aus_idx)
    return p,zin,zau

def robust(p,cap=10):
    p=p.sort_values(["mnn","dist"],ascending=[False,True])
    return p.head(cap)

def zpath(df,muSD):
    x=np.log(df[PATH].clip(lower=1e-9)); return (x-muSD[0])/muSD[1]

def spearman(a,b):
    ra=pd.Series(a).rank().values; rb=pd.Series(b).rank().values
    return float(np.corrcoef(ra,rb)[0,1])

def main():
    global AUS_THR, B, Bc, CI, DR, ELEM, H, HERE, IN_THR, MATCH, NG, NPERM, PATH, PROJ, RES, SEED, _, ac, allp, ar, ar_, c, cidx, da, den, di, el, f, k, m, mang, ng, ngcm_all, ngcm_muSD, ngsa_all, ngsa_muSD, null, ok, ordered, pau_mu, pin_mu, pp, pval, r, ra, rb, rho, rng, rob, sand, use, v, wapp, xi, yi, yiln, zau_ar, zau_pp, zin_ar, zin_pp
    HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ=os.path.dirname(HERE); RES=os.path.join(HERE,"results")
    PROJ=os.path.dirname(HERE); RES=os.path.join(HERE,"results")
    DR=os.path.join(RES,"drainage")
    NPERM = 100000
    SEED = 20260827
    IN_THR=50.0; AUS_THR=25.0   # India >=50% (dense); Australia relaxed to >=25% (sparse, large catchments)
    IN_THR=50.0; AUS_THR=25.0   # India >=50% (dense); Australia relaxed to >=25% (sparse, large catchments)
    CI=dict(La=.237,Yb=.170,Sm=.148,Eu=.0580,Gd=.199)
    MATCH=["Th/Sc","La/Sc","Th/Co","EuEu","La/Yb_n","Nb/Y"]
    PATH=["Zr","Hf","Ti","Ce","Nd","Pr","Dy","P"]
    ar=pd.read_csv(paths.NGCM_TABLE)
    for c in ["Th","Sc","Co","La","Eu","Sm","Gd","Yb","Nb","Y","TiO2","P2O5","Zr","Hf","Ce","Nd","Pr","Dy","LAT","LON"]:
        ar[c]=pd.to_numeric(ar[c],errors="coerce")   # element cols are object (detection-limit strings)
    ar["Ti"]=ar["TiO2"]*1e4*0.5995; ar["P"]=ar["P2O5"]*1e4*0.4364   # oxide wt% -> element mg/kg
    ar["Ti"]=ar["TiO2"]*1e4*0.5995; ar["P"]=ar["P2O5"]*1e4*0.4364   # oxide wt% -> element mg/kg
    ar["k"]=ar.LAT.round(4).astype(str)+"_"+ar.LON.round(4).astype(str)
    sand=india("sandmata"); mang=india("mangalwar")
    sand=india("sandmata"); mang=india("mangalwar")
    NG=paths.NGSA
    with open(NG,encoding="latin-1") as f: H=list(csv.reader(f))[11]
    ELEM=dict(Th=("ICP-MS",),Sc=("ICP-MS",),Nb=("ICP-MS",),Y=("ICP-MS",),La=("ICP-MS",),
              Yb=("ICP-MS",),Ce=("ICP-MS",),Nd=("ICP-MS",),Sm=("ICP-MS",),Eu=("ICP-MS",),
              Gd=("ICP-MS",),Dy=("ICP-MS",),Zr=("ICP-MS",),Hf=("ICP-MS",),Co=("ICP-MS",),
              Pr=("ICP-MS",),Ti=("XRF",),P=("XRF",))
    cidx={el:col(el,m[0]) for el,m in ELEM.items()}
    cidx={k:v for k,v in cidx.items() if v is not None}
    ordered=sorted(cidx.items(),key=lambda kv:kv[1])   # pandas returns usecols in ascending file order
    use=[0]+[v for _,v in ordered]
    ng=pd.read_csv(NG,header=None,skiprows=12,usecols=use,encoding="latin-1",low_memory=False)
    ng.columns=["SITEID"]+[k for k,_ in ordered]
    for c in ng.columns: ng[c]=pd.to_numeric(ng[c],errors="coerce")
    ng=ng.groupby("SITEID").median(numeric_only=True).reset_index()   # one row per site
    wapp=aus(os.path.join(DR,"wa_palaeoprot_contained.csv"),"WA_PP")
    yiln=aus(os.path.join(DR,"yilgarn_contained.csv"),"Y+N")
    ngcm_all=pd.concat([sand,mang]); ngsa_all=pd.concat([wapp,yiln])
    ngcm_all=pd.concat([sand,mang]); ngsa_all=pd.concat([wapp,yiln])
    _,ngcm_muSD=zlog(ngcm_all); _,ngsa_muSD=zlog(ngsa_all)
    _,ngcm_muSD=zlog(ngcm_all); _,ngsa_muSD=zlog(ngsa_all)
    pp,zin_pp,zau_pp=pair(sand,wapp,"Palaeoproterozoic")
    ar_,zin_ar,zau_ar=pair(mang,yiln,"Archaean")
    allp=pd.concat([pp,ar_]).sort_values(["domain","dist"])
    allp.to_csv(os.path.join(RES,"candidate_pairs.csv"),index=False)
    rob=pd.concat([robust(pp),robust(ar_)]).reset_index(drop=True)
    rob.to_csv(os.path.join(RES,"analogue_pairs.csv"),index=False)
    print("candidate pairs: Palaeoprot",len(pp),"(MNN",int(pp.mnn.sum()),") Archaean",len(ar_),"(MNN",int(ar_.mnn.sum()),")")
    print("robust pairs kept:",len(rob),"| median dist",round(rob.dist.median(),2))
    print(rob[["domain","india_sid","aus_sid","dist","mnn"]].to_string(index=False))
    pin_mu=(np.log(ngcm_all[PATH].clip(lower=1e-9)).mean(),np.log(ngcm_all[PATH].clip(lower=1e-9)).std(ddof=0))
    pau_mu=(np.log(ngsa_all[PATH].clip(lower=1e-9)).mean(),np.log(ngsa_all[PATH].clip(lower=1e-9)).std(ddof=0))
    xi,yi=[],[]
    for _,r in rob.iterrows():
        di=ngcm_all[ngcm_all.sid==r.india_sid]; da=ngsa_all[ngsa_all.sid==r.aus_sid]
        if len(di) and len(da):
            xi.append(zpath(di,pin_mu)[PATH].mean(axis=1).values[0])
            yi.append(zpath(da,pau_mu)[PATH].mean(axis=1).values[0])
    xi=np.array(xi); yi=np.array(yi); ok=np.isfinite(xi)&np.isfinite(yi)
    xi=np.array(xi); yi=np.array(yi); ok=np.isfinite(xi)&np.isfinite(yi)
    xi=np.array(xi); yi=np.array(yi); ok=np.isfinite(xi)&np.isfinite(yi)
    xi,yi=xi[ok],yi[ok]
    rho=spearman(xi,yi)
    rng=np.random.default_rng(SEED)
    ra=pd.Series(xi).rank().values; rb=pd.Series(yi).rank().values
    ra=pd.Series(xi).rank().values; rb=pd.Series(yi).rank().values
    B=rb[rng.random((NPERM,len(rb))).argsort(axis=1)]
    Bc=B-B.mean(1,keepdims=True); ac=ra-ra.mean()
    Bc=B-B.mean(1,keepdims=True); ac=ra-ra.mean()
    den=np.sqrt((Bc**2).sum(1)*(ac**2).sum())
    null=np.abs(np.where(den>0,(Bc@ac)/den,np.nan))
    pval=(np.sum(null>=abs(rho))+1)/(NPERM+1)
    print("\nPATHFINDER validation (composite, n=%d pairs): Spearman rho=%.2f, perm p=%.3f"%(len(xi),rho,pval))
    pd.DataFrame({"india":list(rob.india_sid[:len(xi)]),"aus":list(rob.aus_sid[:len(xi)]),
                  "path_india":xi,"path_aus":yi}).to_csv(os.path.join(RES,"pairing_fingerprint_check.csv"),index=False)
    print("wrote analogue pairing csvs")

if __name__ == "__main__":
    main()
