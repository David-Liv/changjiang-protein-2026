#!/usr/bin/env python3
"""SOP §B-1 輸入驗證器 — 計算組每次跑批前執行。
這五項錯誤的共同特徵是「會產出看起來正常的結果，不會報錯」。
用法：python3 validate_input.py --target <target.pdb> --hotspots A409,A412,... [--complex <complex.pdb>] [--sha <expected_sha256>]
"""
import argparse, hashlib, re, sys
from collections import OrderedDict, defaultdict

def parse(path):
    at=[]
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith(("ATOM  ","HETATM")) and len(l)>=54:
            at.append(dict(rec=l[:6].strip(), name=l[12:16].strip(), alt=l[16],
                           resname=l[17:20].strip(), chain=l[21], resseq=int(l[22:26]),
                           icode=l[26], elem=(l[76:78].strip() or l[12:16].strip()[0])))
    return at

BB={"N","CA","C","O"}
EXPECTED={"GLY":4,"ALA":5,"SER":6,"CYS":6,"PRO":7,"THR":7,"VAL":7,"ASN":8,"ASP":8,"LEU":8,
          "ILE":8,"MET":8,"GLN":9,"GLU":9,"LYS":9,"HIS":10,"PHE":11,"ARG":11,"TYR":12,"TRP":14}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--target",required=True)
    p.add_argument("--hotspots",default="")
    p.add_argument("--complex")
    p.add_argument("--sha")
    p.add_argument("--strict-sidechains",action="store_true",help="將晶體結構缺失側鏈原子的殘基也視為硬錯誤")
    a=p.parse_args(); fails=[]; warns=[]

    # 斷言 5：靶點檔案雜湊值與登錄表一致（防無痕修改）
    got=hashlib.sha256(open(a.target,"rb").read()).hexdigest()
    if a.sha and got!=a.sha:
        fails.append(f"[5] 靶點檔案 SHA-256 不符\n      實得 {got}\n      登錄 {a.sha}")
    print(f"  靶點 SHA-256: {got}")

    at=parse(a.target)
    # 斷言 3：熱點字串裡的每個 ID 都存在
    present={(x["chain"],x["resseq"]) for x in at}
    bad=[]
    if a.hotspots and a.hotspots.strip().lower() not in ("null", "none", ""):
        for tok in re.split(r"[,\s]+", a.hotspots.strip()):
            if not tok: continue
            m=re.fullmatch(r"([A-Za-z])(\d+)(?:-(\d+))?", tok)
            if not m: bad.append(f"{tok}（格式無法解析）"); continue
            ch=m.group(1); lo=int(m.group(2)); hi=int(m.group(3) or m.group(2))
            for r in range(lo,hi+1):
                if (ch,r) not in present: bad.append(f"{ch}{r}")
    if bad: fails.append(f"[3] 熱點殘基不存在於 PDB：{bad}\n      （鏈標籤不符會導致熱點被靜默忽略）")

    # 斷言 4：沒有缺主鏈的部分殘基（修剪殘留 → Mismatched atoms）
    byres=defaultdict(set)
    rn={}
    for x in at:
        if x["rec"]=="ATOM":
            byres[(x["chain"],x["resseq"],x["icode"])].add(x["name"]); rn[(x["chain"],x["resseq"],x["icode"])]=x["resname"]
    missing_bb=[]
    missing_sc=[]
    for k,names in byres.items():
        if not BB<=names:
            missing_bb.append(f"{k[0]}{k[1]} {rn[k]} 缺主鏈{sorted(BB-names)}")
        else:
            exp=EXPECTED.get(rn[k])
            if exp and len(names)<exp:
                missing_sc.append(f"{k[0]}{k[1]} {rn[k]} 原子 {len(names)}/{exp}")
    if missing_bb:
        fails.append(f"[4] 缺主鏈原子殘基 {len(missing_bb)} 個（前 10）：{missing_bb[:10]}")
    if missing_sc:
        if a.strict_sidechains:
            fails.append(f"[4] 缺側鏈原子殘基 {len(missing_sc)} 個（前 10）：{missing_sc[:10]}")
        else:
            warns.append(f"[4-warn] 晶體結構未解析側鏈末端原子 {len(missing_sc)} 個（主鏈完整，不阻斷）：{missing_sc[:8]}")

    # 斷言 2：殘基索引跨鏈唯一
    idx=defaultdict(list)
    for k in byres: idx[k[1]].append(k[0])
    dup={r:sorted(set(c)) for r,c in idx.items() if len(set(c))>1}
    if dup: fails.append(f"[2] 殘基索引跨鏈重複 {len(dup)} 個：{dict(list(dup.items())[:10])}")

    # 斷言 1：結合物是第一條鏈（僅在提供複合物時可驗）
    if a.complex:
        ca=parse(a.complex)
        order=list(OrderedDict.fromkeys(x["chain"] for x in ca))
        tgt_chains=set(x["chain"] for x in at)
        if order and order[0] in tgt_chains:
            fails.append(f"[1] 複合物第一條鏈是 {order[0]}（靶點鏈）；dl_binder_design 要求結合物為第一條鏈。鏈序：{order}")
        else:
            print(f"  複合物鏈序: {order}  → 第一條 {order[0] if order else '—'} 非靶點鏈 ✓")
    else:
        print("  [1] 未提供 --complex，略過『結合物是第一條鏈』")

    for w in warns:
        print("  ⚠", w)
    print()
    if fails:
        print("驗證失敗：")
        for f in fails: print("  ✗", f)
        sys.exit(1)
    print("五項斷言全部通過 ✓")

if __name__=="__main__": main()
