import csv, math, random, statistics
from collections import Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                 "figure.dpi":150,"savefig.bbox":"tight"})
INK="#1b2a41"; ACC="#c1502e"; MUT="#8a9aa8"
def sieve(n):
    f=bytearray([1])*(n+1); f[0]=f[1]=0
    i=2
    while i*i<=n:
        if f[i]: f[i*i::i]=bytearray(len(f[i*i::i]))
        i+=1
    return f
C=50000; isp=sieve(C+10)
rows=[(int(r["p"]),int(r["c"])) for r in csv.DictReader(
      open("/mnt/user-data/uploads/e1_N700_C50k__1_.csv"))]
HI=2003
fam={p for p,_ in rows if 100<=p<=HI}
got=Counter(c for p,c in rows if p in fam)
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

fig,(a1,a2)=plt.subplots(1,2,figsize=(7.0,2.7))
a1.hist(sim,bins=44,color=MUT,edgecolor="white",lw=.4)
a1.axvline(mp,color=ACC,lw=1.8)
a1.text(mp-5.3,a1.get_ylim()[1]*0.72,"powers of two\n%.2f  ($p=%.4f$)"%(mp,pval),
        color=ACC,fontsize=8)
a1.set_xlabel("mean primes per offset,\nrandom 15-offset subsets")
a1.set_ylabel("count")
a1.set_title("(a) the original comparison",fontsize=9,loc="left")

# replication: matched pairs
runs=[("$p\\leq157$",1864,146,1857,143),("$p\\leq419$",16305,634,16285,605),
      ("$p\\leq797$",47323,936,47301,927)]
labs=[r[0] for r in runs]+["pooled"]
ac=sum(r[1] for r in runs); ah=sum(r[2] for r in runs)
bc=sum(r[3] for r in runs); bh=sum(r[4] for r in runs)
rat=[];err=[]
for _,c1,h1,c2,h2 in runs+[("pooled",ac,ah,bc,bh)]:
    rat.append((h1/c1)/(h2/c2)); err.append(math.sqrt(1/h1+1/h2))
y=range(len(labs))
a2.axvline(1.0,color=INK,lw=1,zorder=1)
a2.axvline(1.23,color=ACC,lw=1.4,ls="--",zorder=1)
a2.errorbar(rat,y,xerr=[1.96*e*r for e,r in zip(err,rat)],fmt="o",color=INK,
            ms=4,lw=1.2,capsize=3,zorder=3)
a2.set_yticks(list(y)); a2.set_yticklabels(labs,fontsize=8)
a2.set_xlim(0.80,1.35); a2.invert_yaxis()
a2.text(1.20,2.6,"claimed\n+23%",color=ACC,fontsize=8,ha="left")
a2.set_xlabel("primes with $\\Pr=1$ / primes with $\\Pr$ prime\n(matched in family, size and count)")
a2.set_title("(b) the matched replication",fontsize=9,loc="left")
fig.savefig("fig/perm.pdf"); fig.savefig("fig/perm.eps")
print("ok  pooled ratio %.4f"%(rat[-1]))
