#!/usr/bin/env python3
"""SOP §C-2 組成檢查。輸入序列（FASTA 或純序列，每行一條），輸出判讀欄位。
門檻說明：SOP 的 E／K 比例 ROC AUC（0.77／0.73）來自 Round 2 的 378 條實測，
公開資料沒有給數值門檻，只給方向（高者佳）。故本腳本輸出數值與相對分位，不做硬性通過／退回。
pH 6.0 註記：本靶點的親和力在 pH 6.0 測定，故同時輸出該 pH 下的淨電荷。"""
import sys, re
PKA_POS={"K":10.54,"R":12.48,"H":6.04}; PKA_NEG={"D":3.90,"E":4.07,"Y":10.46,"C":8.18}
def charge(s,pH,freeC=True):
    q=1/(1+10**(pH-8.60))-1/(1+10**(3.65-pH))
    for aa,pk in PKA_POS.items(): q+=s.count(aa)/(1+10**(pH-pk))
    for aa,pk in PKA_NEG.items():
        if aa=="C" and not freeC: continue
        q-=s.count(aa)/(1+10**(pk-pH))
    return q
MW={'A':71.08,'R':156.19,'N':114.10,'D':115.09,'C':103.14,'Q':128.13,'E':129.12,'G':57.05,
    'H':137.14,'I':113.16,'L':113.16,'K':128.17,'M':131.19,'F':147.18,'P':97.12,'S':87.08,
    'T':101.10,'V':99.13,'W':186.21,'Y':163.18}
HYD=set("AVLIMFWYC")
def report(name,s):
    s=s.upper(); n=len(s)
    e,k=s.count("E")/n, s.count("K")/n
    cys=s.count("C")
    runs=max((len(m.group(0)) for m in re.finditer(r"[AVLIMFWYC]+",s)), default=0)
    homo=max((len(m.group(0)) for m in re.finditer(r"(.)\1+",s)), default=1)
    mw=sum(MW.get(c,110) for c in s)/1000+0.018
    print(f"\n{name}  長度 {n}  估算分子量 {mw:.1f} kDa")
    print(f"  E 比例 {e*100:5.1f}%   K 比例 {k*100:5.1f}%   E+K {(e+k)*100:5.1f}%      （ROC AUC 0.77／0.73，高者佳）")
    print(f"  淨電荷  pH7.4 {charge(s,7.4,cys%2==1):+6.1f}   pH6.0 {charge(s,6.0,cys%2==1):+6.1f}   His {s.count('H')}")
    print(f"  Cys {cys} → {'偶數（可配對）' if cys%2==0 else '★奇數：必有落單 Cys，退回重做該位置'}")
    print(f"  最長連續疏水片段 {runs}   疏水殘基總比 {100*sum(1 for c in s if c in HYD)/n:.1f}%   最長同聚物串 {homo}")
    flags=[]
    if cys%2: flags.append("Cys 奇數 → 退回（SOP §C-5）")
    if n>250: flags.append("長度 >250，超出提交上限")
    if n<10: flags.append("長度 <10，低於提交下限")
    if runs>=8: flags.append(f"連續疏水 {runs} 殘基，聚集風險")
    if homo>=5: flags.append(f"同聚物串 {homo}，表達與合成風險")
    print(f"  旗標：{flags or '無'}")
def main():
    txt=sys.stdin.read() if len(sys.argv)<2 else open(sys.argv[1]).read()
    name=None; buf=[]
    for line in txt.splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if buf: report(name or "seq","".join(buf)); buf=[]
            name=line[1:]
        else: buf.append(line)
    if buf: report(name or "seq","".join(buf))
if __name__=="__main__": main()
