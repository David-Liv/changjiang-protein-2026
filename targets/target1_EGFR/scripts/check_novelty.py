"""SOP §C-1 新穎性硬閘檢查器。
Needleman–Wunsch 全域比對（純標準庫，BLOSUM62 簡化為 match/mismatch + gap）。
輸出：全長同一性、框架區／CDR 分區同一性、相差胺基酸數。
"""
import sys, re
def read_fasta(p):
    out={}; k=None; b=[]
    for l in open(p):
        l=l.strip()
        if not l: continue
        if l[0]==">":
            if k: out[k]="".join(b)
            k=l[1:].split()[0]; b=[]
        else: b.append(l)
    if k: out[k]="".join(b)
    return out
def nw(a,b,match=2,mis=-1,gap=-2):
    n,m=len(a),len(b)
    F=[[0]*(m+1) for _ in range(n+1)]; P=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): F[i][0]=gap*i; P[i][0]=1
    for j in range(1,m+1): F[0][j]=gap*j; P[0][j]=2
    for i in range(1,n+1):
        ai=a[i-1]
        for j in range(1,m+1):
            d=F[i-1][j-1]+(match if ai==b[j-1] else mis)
            u=F[i-1][j]+gap; l=F[i][j-1]+gap
            best=max(d,u,l); F[i][j]=best
            P[i][j]=0 if best==d else (1 if best==u else 2)
    i,j=n,m; A=[];B=[]
    while i>0 or j>0:
        p=P[i][j]
        if p==0: A.append(a[i-1]); B.append(b[j-1]); i-=1; j-=1
        elif p==1: A.append(a[i-1]); B.append("-"); i-=1
        else: A.append("-"); B.append(b[j-1]); j-=1
    return "".join(reversed(A)),"".join(reversed(B))
def ident(a,b):
    A,B=nw(a,b)
    aln=sum(1 for x,y in zip(A,B) if x!="-" and y!="-")
    same=sum(1 for x,y in zip(A,B) if x==y and x!="-")
    diff=sum(1 for x,y in zip(A,B) if x!=y)          # 含 gap 的位置
    return same, aln, 100*same/max(aln,1), diff

def cdr3(s):
    """VHH CDR3：最後一個 C 與 WG.G 之間。"""
    m=re.search(r"C(?!.*C)(.*?)(WG.G)", s, re.S)
    return m.group(1) if m else None

def main(path):
    seqs=read_fasta(path)
    vhh={k:v for k,v in seqs.items() if "VHH" in k.upper()}
    print("== 排除集內容 ==")
    CANON="QVQLVESGGGLVQPGGSLRLSCAAS"
    for k,v in seqs.items():
        c3=cdr3(v)
        print(f"  {k:<18} 長度 {len(v):>3}  Cys {v.count('C')}  CDR3 {'('+str(len(c3))+' aa) '+c3 if c3 else '—'}")
    print("\n== 缺口偵測：比對 VHH 典型框架 FR1 開頭 ==")
    print(f"  典型 FR1 : {CANON}")
    for k,v in vhh.items():
        print(f"  {k:<12}: {v[:25]}   {'★ 明顯短缺' if len(v[:25].replace('-',''))<25 or 'SCAAS' in v[:16] else ''}")
    print("\n== 已知 VHH 兩兩同一性（這決定新穎性閘有多緊）==")
    ks=list(vhh)
    print(f"  {'':<14}"+"".join(f"{k:>14}" for k in ks))
    for a in ks:
        row=f"  {a:<14}"
        for b in ks:
            if a==b: row+=f"{'—':>14}"
            else:
                s_,al,pct,_=ident(vhh[a],vhh[b]); row+=f"{pct:>13.1f}%"
        print(row)
    print("\n== 對框架區與 CDR3 分開看 ==")
    for a in ks:
        for b in ks:
            if a>=b: continue
            c3a,c3b=cdr3(vhh[a]),cdr3(vhh[b])
            _,_,full,_=ident(vhh[a],vhh[b])
            if c3a and c3b:
                _,_,c3pct,_=ident(c3a,c3b)
                fa=vhh[a].replace(c3a,""); fb=vhh[b].replace(c3b,"")
                _,_,fwpct,_=ident(fa,fb)
                print(f"  {a} vs {b}:  全長 {full:.1f}%   框架 {fwpct:.1f}%   CDR3 {c3pct:.1f}%")
if __name__=="__main__": main(sys.argv[1] if len(sys.argv)>1 else "data/raw/known_binders.fasta")
