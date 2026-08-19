import csv, math, random, statistics
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                 "figure.dpi":150,"savefig.bbox":"tight"})
EG=math.exp(0.5772156649)
INK="#1b2a41"; ACC="#c1502e"; MUT="#8a9aa8"

def sieve(n):
    f=bytearray([1])*(n+1); f[0]=f[1]=0
    i=2
    while i*i<=n:
        if f[i]: f[i*i::i]=bytearray(len(f[i*i::i]))
        i+=1
    return f

C=50000
isp=sieve(C+10)
src="/mnt/user-data/uploads/e1_N700_C50k__1_.csv"
rows=[(int(r["p"]),int(r["digits"]),int(r["c"])) for r in csv.DictReader(open(src))]
cnt=Counter(p for p,_,_ in rows); dig={p:d for p,d,_ in rows}
prs=[q for q in range(3,5300) if isp[q]]

# ---------- Fig 1: admissible offsets, and how the structure thins ----------
adm={}
for p in prs:
    n=0
    for c in range(2,C+1,2):
        m=c
        while m%2==0: m//=2
        if m==1 or (isp[m] and m>p): n+=1
    adm[p]=n
fig,ax=plt.subplots(figsize=(5.0,2.9))
xs=sorted(adm); ax.plot(xs,[adm[p] for p in xs],color=INK,lw=1.4)
ax.set_xscale("log"); ax.set_xlabel("$p$"); ax.set_ylabel("admissible offsets $c\\leq 50000$")
ax.axhline(math.log2(C),color=ACC,ls="--",lw=1)
ax.text(3.5,math.log2(C)+120,"the %d pure powers of two survive every $p$"%int(math.log2(C)),
        color=ACC,fontsize=7.5)
ax.set_title("Admissibility thins the offset axis as $p$ grows",fontsize=9.5,loc="left")
fig.savefig("fig/admissible.pdf")

# ---------- Fig 2: Bateman-Horn check ----------
ratios=[]
for p in prs:
    if p not in dig: continue
    exp=2*adm[p]*EG*math.log(p)/(dig[p]*math.log(10))
    if exp>5: ratios.append(cnt[p]/exp)
fig,ax=plt.subplots(figsize=(5.0,2.9))
ax.hist(ratios,bins=34,color=MUT,edgecolor="white",lw=.5)
m=statistics.mean(ratios)
ax.axvline(1.0,color=INK,lw=1.2,label="prediction")
ax.axvline(m,color=ACC,lw=1.6,ls="--",label="observed mean %.3f"%m)
ax.set_xlabel("observed primes / predicted primes"); ax.set_ylabel("families")
ax.legend(frameon=False,fontsize=8)
ax.set_title("Bateman--Horn against %d families, %s primes"%(len(ratios),format(len(rows),",")),
             fontsize=9.5,loc="left")
fig.savefig("fig/batemanhorn.pdf")

# ---------- Fig 3: permutation test for the power-of-two effect ----------
HI=2003
fam={p for p,_,_ in rows if 100<=p<=HI}
got=Counter(c for p,_,c in rows if p in fam)
pw=[];oth=[]
for c in range(2,C+1,2):
    m=c
    while m%2==0: m//=2
    if m==1: pw.append(got.get(c,0))
    elif isp[m] and m>HI: oth.append(got.get(c,0))
mp=statistics.mean(pw)
random.seed(1)
sim=[statistics.mean(random.sample(oth,len(pw))) for _ in range(20000)]
pval=sum(1 for s in sim if s>=mp)/len(sim)
fig,ax=plt.subplots(figsize=(5.0,2.9))
ax.hist(sim,bins=48,color=MUT,edgecolor="white",lw=.4)
ax.axvline(mp,color=ACC,lw=1.8)
ax.annotate("powers of two\n%.2f  ($p=%.4f$)"%(mp,pval),xy=(mp,0),
            xytext=(mp-4.6,ax.get_ylim()[1]*0.80),color=ACC,fontsize=8,
            arrowprops=dict(arrowstyle="->",color=ACC,lw=1))
ax.set_xlabel("mean primes per offset, random 15-offset subsets")
ax.set_ylabel("count")
ax.set_title("The one effect that does not average away",fontsize=9.5,loc="left")
fig.savefig("fig/perm.pdf")

# ---------- Fig 4: Poisson check on the two large runs ----------
from math import exp as e_, factorial
fig,axes=plt.subplots(1,2,figsize=(6.4,2.6))
for ax,(lab,lam,obs) in zip(axes,[("6003 digits",1.70,4),("20501 digits",2.81,1)]):
    ks=list(range(0,11))
    ax.bar(ks,[e_(-lam)*lam**k/factorial(k) for k in ks],color=MUT,edgecolor="white",lw=.5)
    ax.bar([obs],[e_(-lam)*lam**obs/factorial(obs)],color=ACC,edgecolor="white",lw=.5)
    ax.set_title("%s   $\\lambda=%.2f$, found %d"%(lab,lam,obs),fontsize=8.5,loc="left")
    ax.set_xlabel("primes found"); ax.set_xticks(range(0,11,2))
axes[0].set_ylabel("probability")
fig.savefig("fig/poisson.pdf")

# ---------- Fig 5: cost model ----------
D=[2000,3000,5000,8000,12000,16000,20500]
a=2.5e-9
fig,ax=plt.subplots(figsize=(5.0,2.9))
ax.loglog(D,[a*d**2.24 for d in D],color=INK,lw=1.5,label="$1.86\\times10^{-9}\\,D^{2.24}$")
ax.loglog([2000,12002,20501],[0.060,3.42,11.4],"o",color=ACC,ms=5,label="measured")
ax.set_xlabel("decimal digits $D$"); ax.set_ylabel("seconds per test, one core")
ax.legend(frameon=False,fontsize=8)
ax.set_title("Cost of one strong probable-prime test  ($R^2=0.983$)",fontsize=9.5,loc="left")
fig.savefig("fig/cost.pdf")
print("figures written; permutation p =",pval,"; BH mean =",round(m,3))
