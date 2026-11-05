#!/usr/bin/env python3
"""TNF-α ipSAE 測試管線 —— 單檔版（只依賴 numpy；摺疊才需要 torch）。

子命令:
  selftest  不需 GPU 的自測：ipSAE 性質、PDB 分鏈改寫、介面指標、triage
  triage    序列預檢（長度 ≤ 250 aa、Cysteine 奇偶、淨電荷、GRAVY、疏水連續段）
  fold      ESMFold 推理（HuggingFace 後端，不需編譯 openfold）
  fold-esm  ESMFold 推理（fair-esm 後端，ModelScope / 24GB 卡路線）
  score     把輸出全部過 ipSAE + 介面指標 + triage，匯總成兩張表與 CSV
  compare   與晶體結構疊合（自動從 RCSB 下載、自動對位）

用法:
  python tnf_pipeline.py selftest
  python tnf_pipeline.py triage --spec spec_stage23.tsv
  python tnf_pipeline.py fold  --spec spec_stage1_mono.tsv --out runs/s1 --chunk 64
  python tnf_pipeline.py score --dir runs/s1 --cutoffs 10 15 --csv runs/s1/summary.csv
  python tnf_pipeline.py compare runs/s1/1a_TNF_mono__r4L25.pdb --fetch 5M2J

spec 檔每行: 名稱<TAB>鏈1序列:鏈2序列[:鏈3序列...]（第一條鏈視為 binder）
"""
import argparse, csv, glob, io, json, math, os, sys, time, urllib.request
import numpy as np

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
JUMP = 512          # 鏈交界處的 position_ids 跳號（社群標準作法）


# ═══════════════════════════════════════════════════════════════════
#  1. ipSAE 計分（Dunbrack 定義）
# ═══════════════════════════════════════════════════════════════════
# 對每個配對 (i∈鏈A, j∈鏈B) 取其 PAE 值 d，僅計入 PAE < cutoff 的 j，令 n = 通過者數。
#   d0 = max(1.24 * (L - 15)^(1/3) - 1.8, 1.0)，其中 L = max(n, 27)
#   （原版 calc_d0 在 L < 27 時夾到 27；與單純夾 1.0 Å 差 < 4%，此處照原版寫法。）
#   每個 i 的分數 = 對通過的 j 取 1/(1+(d/d0)^2) 的平均；ipSAE(A→B) = 對所有 i 取最大值。
# 非對稱：ipSAE_max = max(A→B, B→A)；ipSAE_min = min(...)（binder design 慣例用 min）。
# 注意：不同 PAE cutoff 下的數值不可互相比較。官方 cutoff 未公開，預設 10 Å。

def d0_from_n(n):
    """純量進→float 出；陣列進→陣列出。"""
    arr = np.asarray(n, dtype=float)
    d0 = 1.24 * np.cbrt(np.maximum(arr, 27.0) - 15.0) - 1.8
    d0 = np.maximum(d0, 1.0)
    return float(d0) if np.ndim(n) == 0 else d0


def ipsae_directional(pae, idx_from, idx_to, cutoff):
    """ipSAE(from → to)。pae[i, j] = 以 i 為對位殘基時 j 的預測對位誤差。

    NumPy 向量化：一次算完所有 i，批次篩選上千條序列時比逐列迴圈快一個量級。
    回傳 (best_score, best_i 全域索引, best_n, best_partners 全域索引)。
    """
    idx_from = np.asarray(idx_from, dtype=int)
    idx_to = np.asarray(idx_to, dtype=int)
    if idx_from.size == 0 or idx_to.size == 0:
        return 0.0, None, 0, np.array([], dtype=int)
    sub = pae[np.ix_(idx_from, idx_to)].astype(float)          # (F, T)
    keep = sub < cutoff
    n = keep.sum(axis=1)
    d0 = d0_from_n(n)
    w = np.where(keep, 1.0 / (1.0 + (sub / d0[:, None]) ** 2), 0.0)
    score = np.where(n > 0, w.sum(axis=1) / np.maximum(n, 1), 0.0)
    b = int(np.argmax(score))
    if score[b] <= 0.0:
        return 0.0, None, 0, np.array([], dtype=int)
    return float(score[b]), int(idx_from[b]), int(n[b]), idx_to[keep[b]]


def _entry(pae, ia, ib, cutoff):
    ab, i_a, na, _ = ipsae_directional(pae, ia, ib, cutoff)
    ba, i_b, nb, _ = ipsae_directional(pae, ib, ia, cutoff)
    return {"asym_a_to_b": round(ab, 4), "asym_b_to_a": round(ba, 4),
            "max": round(max(ab, ba), 4), "min": round(min(ab, ba), 4),
            "n_pass_a": na, "n_pass_b": nb,
            "best_res_a": i_a, "best_res_b": i_b, "cutoff": cutoff}


def ipsae(pae, chain_of, cutoff=10.0, binder=None, merge_target=True):
    """每個鏈對的 asym / max / min；鏈數 > 2 時額外給 binder 對合併靶點的讀數。

    多聚體靶點（binder A + TNF 三聚體 B/C/D）：pairwise 會把跨亞基界面切成兩半，
    使每一半的 n 變小、d0 懲罰過重而低估分數。"A-ALL" 把 B+C+D 視為單一靶點整體
    計算 n 與 d0 —— 這已偏離原版定義，僅作診斷用，不可當成官方分數。
    """
    lab = np.asarray(chain_of)
    chains = sorted(set(lab.tolist()))
    idx = {c: np.where(lab == c)[0] for c in chains}
    out = {}
    for a in chains:
        for b in chains:
            if a < b:
                out[f"{a}-{b}"] = _entry(pae, idx[a], idx[b], cutoff)
    if len(chains) > 2 and merge_target:
        bc = binder if binder in chains else chains[0]
        others = [c for c in chains if c != bc]
        e = _entry(pae, idx[bc], np.concatenate([idx[c] for c in others]), cutoff)
        e["merged_chains"] = "".join(others)
        out[f"{bc}-ALL"] = e
    return out


def summarize(res, binder="A"):
    """挑出 binder–target 的 pair_max（所有含 binder 的鏈對最大值）與合併讀數。
    這是修掉「只看字母序第一個鏈對 A-B」的地方：跨亞基結合可能落在 A-C 或 A-D。"""
    pairs = {k: v for k, v in res.items()
             if "-ALL" not in k and binder in k.split("-")}
    if not pairs:
        pairs = {k: v for k, v in res.items() if "-ALL" not in k}
    if not pairs:
        return None
    key = max(pairs, key=lambda k: pairs[k]["max"])
    merged = res.get(f"{binder}-ALL")
    return {"pair": key, "pair_max": pairs[key]["max"], "pair_min": pairs[key]["min"],
            "n_pass": max(pairs[key]["n_pass_a"], pairs[key]["n_pass_b"]),
            "all_pairs": {k: v["max"] for k, v in sorted(pairs.items())},
            "merged_max": merged["max"] if merged else None,
            "merged_min": merged["min"] if merged else None,
            "merged_n": max(merged["n_pass_a"], merged["n_pass_b"]) if merged else None}


# ═══════════════════════════════════════════════════════════════════
#  2. PDB 讀寫與清理
# ═══════════════════════════════════════════════════════════════════

def read_pdb(path):
    """回傳 {chain: [residue, ...]}，residue = dict(resseq, aa, atoms, rep, ca, bfac)。
    rep = CB 座標（Gly 用 CA）；atoms 為該殘基所有重原子座標。"""
    chains, cur = {}, {}
    for line in open(path, encoding="utf-8"):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        res3 = line[17:20].strip()
        if res3 not in AA3 or line[16] not in (" ", "A"):
            continue
        elem = line[76:78].strip() or line[12:16].strip()[:1]
        if elem == "H":
            continue
        ch, key = line[21], (line[21], line[22:27])
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        try:
            b = float(line[60:66])
        except ValueError:
            b = 0.0
        r = cur.get(key)
        if r is None:
            r = cur[key] = {"resseq": int(line[22:26]), "aa": AA3[res3],
                            "atoms": [], "rep": None, "ca": None, "b": []}
            chains.setdefault(ch, []).append(r)
        r["atoms"].append(xyz); r["b"].append(b)
        atom = line[12:16].strip()
        if atom == "CB":
            r["rep"] = xyz
        elif atom == "CA":
            r["ca"] = xyz
    for rs in chains.values():
        for r in rs:
            if r["rep"] is None:
                r["rep"] = r["ca"]
            r["atoms"] = np.array(r["atoms"])
            r["bfac"] = float(np.mean(r["b"])) if r["b"] else 0.0
    return chains


def clean_pdb(raw_pdb, lab):
    """剔除 linker 原子，並把鏈 ID 與殘基編號改寫成 lab 指定的鏈、各鏈從 1 起算。

    為什麼不讀 resSeq 欄：HF 的 output_to_pdb 把帶 +512 跳號的 position_ids 直接
    寫進該欄（四鏈時會長成 1..115, 116..140, 653..809 …），拿它當索引會錯位。
    改依原子記錄的出現順序切殘基群組，第 k 群對應 lab[k]，與編號數值無關。
    """
    order, seen = [], {}
    for line in raw_pdb.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21], line[22:27])
            if key not in seen:
                seen[key] = len(order); order.append(key)
    if len(order) != len(lab):
        raise ValueError(f"PDB 殘基數 {len(order)} 與鏈標籤數 {len(lab)} 不符，"
                         f"無法安全改寫編號")
    newmap, counters = {}, {}
    for k, key in enumerate(order):
        ch = lab[k]
        if ch is None:
            continue                       # linker：整個殘基丟掉
        counters[ch] = counters.get(ch, 0) + 1
        newmap[key] = (ch, counters[ch])

    out, serial, prev = [], 0, None
    for line in raw_pdb.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        nm = newmap.get((line[21], line[22:27]))
        if nm is None:
            continue
        ch, rs = nm
        resname = line[17:20]
        if prev is not None and ch != prev[0]:
            serial += 1
            out.append(f"TER   {serial:5d}      {prev[2]:>3s} {prev[0]}{prev[1]:4d}")
        serial += 1
        out.append(f"{line[:6]}{serial:5d}{line[11:21]}{ch}{rs:4d} {line[27:]}")
        prev = (ch, rs, resname)
    if prev is not None:
        serial += 1
        out.append(f"TER   {serial:5d}      {prev[2]:>3s} {prev[0]}{prev[1]:4d}")
    out.append("END")
    return "\n".join(out) + "\n"


def plddt_per_res(pl):
    """壓成每殘基一個值並統一成 0–100 尺度。
    ESMFold 的 plddt 是 (L, 37) 每原子值：取 CA（atom37 索引 1）而非對 37 個原子
    取平均，避免把不存在的原子算進去。"""
    pl = np.asarray(pl, dtype=float)
    if pl.ndim > 1:
        pl = pl[..., 1] if pl.shape[-1] in (14, 37) else pl.mean(-1)
    if pl.size and np.nanmax(pl) <= 1.5:
        pl = pl * 100.0
    return pl


def build(pairs, linker_len=25):
    """pairs = [(鏈ID, 序列), ...]（依實際輸入模型的順序）。
    回傳 (串接序列, position_ids, 每個位置的鏈 ID；linker 位置為 None)。"""
    linker = "G" * linker_len
    full, pos, lab, p = "", [], [], 0
    for k, (cid, s) in enumerate(pairs):
        if k > 0:
            full += linker
            for _ in linker:
                pos.append(p); lab.append(None); p += 1
            p += JUMP
        for _ in s:
            pos.append(p); lab.append(cid); p += 1
        full += s
    return full, np.array(pos), lab


# ═══════════════════════════════════════════════════════════════════
#  3. 介面幾何指標與 pDockQ
# ═══════════════════════════════════════════════════════════════════
# ipSAE 只讀 PAE，PAE 低不代表兩條鏈真的靠在一起。陽性對照全部 < 0.2 時，
# 接觸數可一眼分辨「cutoff 太嚴格」與「根本沒對接上」。
# pDockQ 依 Bryant, Pozzati & Elofsson (2022) Nat Commun 13:1265：
#   n_contacts = CB–CB（Gly 用 CA）< 8 Å 的跨鏈殘基對數；界面殘基 = 參與接觸者
#   x = mean(界面殘基 pLDDT) * log(n_contacts)
#   pDockQ = 0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018
#   判讀：> 0.23 約為可接受模型，> 0.5 多為正確界面。
# pLDDT 由 PDB 的 B-factor 欄讀入（ESMFold 把 pLDDT 寫在該欄）。

def _stack(chains, ids):
    reps, bs, atoms = [], [], []
    for c in ids:
        for r in chains.get(c, []):
            if r["rep"] is None:
                continue
            reps.append(r["rep"]); bs.append(r["bfac"]); atoms.append(r["atoms"])
    return (np.array(reps) if reps else np.zeros((0, 3))), np.array(bs), atoms


def interface(pdb_path, binder="A", target=None, contact_cutoff=8.0):
    chains = read_pdb(pdb_path)
    if binder not in chains:
        return None
    tgt = list(target) if target else [c for c in sorted(chains) if c != binder]
    if not tgt:
        return None
    A, bA, atomsA = _stack(chains, [binder])
    B, bB, atomsB = _stack(chains, tgt)
    if len(A) == 0 or len(B) == 0:
        return None

    d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    hit = d < contact_cutoff
    n_contacts = int(hit.sum())
    ia, ib = np.where(hit.any(axis=1))[0], np.where(hit.any(axis=0))[0]

    min_heavy = float("inf")           # 最小重原子距離：不受 CB 近似影響
    for i in ia:
        for j in ib:
            dd = np.linalg.norm(atomsA[i][:, None, :] - atomsB[j][None, :, :],
                                axis=-1).min()
            min_heavy = min(min_heavy, float(dd))
    if not np.isfinite(min_heavy):
        min_heavy = float(d.min())

    pl = np.concatenate([bA[ia], bB[ib]]) if n_contacts else np.array([])
    if pl.size and np.nanmax(pl) <= 1.5:
        pl = pl * 100.0
    if n_contacts > 0 and pl.size:
        x = float(pl.mean()) * math.log(n_contacts)
        pdockq = 0.724 / (1.0 + math.exp(-0.052 * (x - 152.611))) + 0.018
    else:
        pdockq = 0.0
    return {"n_contacts_cb8": n_contacts, "n_iface_res": int(ia.size + ib.size),
            "min_cb_dist": round(float(d.min()), 2),
            "min_heavy_dist": round(min_heavy, 2),
            "iface_plddt": round(float(pl.mean()), 1) if pl.size else 0.0,
            "pdockq": round(pdockq, 4), "target_chains": "".join(tgt)}


# ═══════════════════════════════════════════════════════════════════
#  4. 序列 triage（大賽硬性約束）
# ═══════════════════════════════════════════════════════════════════
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


def triage(seq, max_len=250):
    """- 長度 ≤ max_len（預設 250 aa）
    - Cysteine 計數；奇數個 C 代表有未配對半胱氨酸 → 錯誤摺疊／二聚體風險
    - 淨電荷（pH 7 近似：K+R+0.1*H 減 D+E）
    - GRAVY（Kyte–Doolittle）；> 0.1 偏疏水 → 聚集風險（奈米抗體常態 -0.4 ~ -0.1）
    - 最長連續疏水片段"""
    s = (seq or "").upper()
    n = len(s)
    nc = s.count("C")
    pos = s.count("K") + s.count("R") + 0.1 * s.count("H")
    neg = s.count("D") + s.count("E")
    gravy = sum(KD.get(c, 0.0) for c in s) / n if n else 0.0
    run = best = 0
    for c in s:
        run = run + 1 if KD.get(c, -9) > 1.5 else 0
        best = max(best, run)
    flags = []
    if n > max_len:
        flags.append(f"長度{n}>{max_len}")
    if nc % 2 == 1:
        flags.append(f"奇數C({nc})")
    if gravy > 0.1:
        flags.append(f"GRAVY{gravy:+.2f}")
    if best >= 6:
        flags.append(f"疏水連續{best}")
    if any(c not in KD for c in s):
        flags.append("非標準殘基")
    return {"len": n, "n_cys": nc, "cys_odd": nc % 2 == 1,
            "net_charge": round(pos - neg, 1), "gravy": round(gravy, 3),
            "max_hydrophobic_run": best, "flags": flags, "pass": not flags}


def read_spec(path):
    """回傳 [(名稱, [(鏈ID, 序列), ...]), ...]；第一條鏈為 binder。"""
    jobs = []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, spec = line.split("\t")
        seqs = spec.split(":")
        jobs.append((name, [(chr(ord("A") + k), s) for k, s in enumerate(seqs)]))
    return jobs


def cmd_triage(a):
    specs = a.spec or sorted(glob.glob("spec_stage*.tsv"))
    if a.seq:
        print(json.dumps(triage(a.seq, a.max_len), indent=1, ensure_ascii=False)); return
    if not specs:
        sys.exit("找不到 spec 檔；用 --spec 指定或 --seq 直接給序列")
    print(f"{'spec':<26}{'名稱':<30}{'長度':>5}{'C數':>4}{'淨電荷':>8}{'GRAVY':>8}  旗標")
    bad = 0
    for f in specs:
        for name, pairs in read_spec(f):
            t = triage(pairs[0][1], a.max_len)        # 第一條鏈 = binder
            flags = "、".join(t["flags"]) if t["flags"] else "OK"
            bad += 0 if t["pass"] else 1
            print(f"{f:<26}{name:<30}{t['len']:>5}{t['n_cys']:>4}"
                  f"{t['net_charge']:>+8.1f}{t['gravy']:>+8.2f}  {flags}")
    print(f"\n{bad} 條序列有旗標；長度超標或奇數 C 的序列不該浪費 GPU，先修再摺疊。")


# ═══════════════════════════════════════════════════════════════════
#  5. 摺疊：HuggingFace 後端
# ═══════════════════════════════════════════════════════════════════
# 為什麼沒有 --seeds：ESMFold 在 model.eval() 下沒有 dropout、也沒有 MSA 隨機採樣，
# 整條推理是確定性的 —— 同序列跑不同 torch.manual_seed 會得到完全相同的 PAE
# （標準差 0，假穩定）。要測穩定度請改變真正影響推理的東西：
#   --recycles 1 2 3 4 / --linkers 20 25 30 / --chain-orders fwd rev
# --check-determinism 可當場驗證「同設定跑兩次 PAE 最大差 ≈ 0」。

def _seed_warning(a):
    import torch
    if a.seeds and len(a.seeds) > 1:
        print(f"⚠ --seeds {a.seeds} 被忽略為單次執行：ESMFold 在 eval() 下沒有 dropout "
              f"與 MSA 採樣，不同 seed 會給出完全相同的 PAE。\n"
              f"  要測穩定度請用 --recycles / --linkers / --chain-orders。")
    seed = (a.seeds or [1])[0]
    torch.manual_seed(seed); np.random.seed(seed)


def pick_device(req):
    import torch
    if req != "auto":
        return req
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


_HF_CACHE = {}


def load_hf(device, chunk, fp16):
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    key = (device, bool(fp16))
    if key in _HF_CACHE:
        tok, m = _HF_CACHE[key]
        if device == "cuda":
            m.trunk.set_chunk_size(chunk or None)
        return tok, m
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    m = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    m = m.eval()
    if device == "cuda":
        m = m.cuda()
        if fp16:
            m.esm = m.esm.half()
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.backends.cudnn.benchmark = True
    elif device == "mps":
        m = m.to("mps")          # Mac GPU：保持 fp32，fp16 在 MPS 上數值不穩
    m.trunk.set_chunk_size(chunk or None)
    _HF_CACHE[key] = (tok, m)
    return tok, m


def infer_hf(model, tok, full, pos, device, recycles):
    import torch
    enc = tok([full], return_tensors="pt", add_special_tokens=False)
    enc = {k: v.to(device) for k, v in enc.items()}
    enc["position_ids"] = torch.tensor(pos, dtype=torch.long, device=device).unsqueeze(0)
    with torch.no_grad():
        out = model(**enc, num_recycles=recycles)
    return (out.predicted_aligned_error[0].float().cpu().numpy(),
            plddt_per_res(out.plddt[0].float().cpu().numpy()),
            model.output_to_pdb(out)[0])


def infer_hf_batch(model, tok, full_list, pos_list, device, recycles):
    """對相同總長度的序列進行 Batch 並行推理（無 Padding 零誤差，壓榨 A100 30–36GB 顯存）。"""
    import torch
    if len(full_list) == 1:
        pae, pl, raw = infer_hf(model, tok, full_list[0], pos_list[0], device, recycles)
        return [(pae, pl, raw)]
    enc = tok(full_list, return_tensors="pt", add_special_tokens=False)
    enc = {k: v.to(device) for k, v in enc.items()}
    enc["position_ids"] = torch.tensor(np.stack(pos_list, axis=0), dtype=torch.long, device=device)
    with torch.no_grad():
        out = model(**enc, num_recycles=recycles)
    pdbs = model.output_to_pdb(out)
    results = []
    for b in range(len(full_list)):
        pae_b = out.predicted_aligned_error[b].float().cpu().numpy()
        pl_b = plddt_per_res(out.plddt[b].float().cpu().numpy())
        results.append((pae_b, pl_b, pdbs[b]))
    return results


def _auto_batch_size(device, seq_len):
    """根據 A100 顯存自動計算並行 Batch 大小，將顯存利用率提升至 30–36 GB。"""
    if device != "cuda":
        return 1
    try:
        import torch
        gb = torch.cuda.get_device_properties(0).total_memory / 1e9
        if gb >= 70:
            return 8 if seq_len <= 320 else (2 if seq_len <= 700 else 1)
        if gb >= 38:
            return 4 if seq_len <= 320 else (1 if seq_len > 450 else 2)
    except Exception:
        pass
    return 1


def _is_tag_completed(outdir, tag):
    """檢查單條序列變體是否已經完整落盤（支援斷點續跑）。"""
    for ext in (".pdb", ".pae.npy", ".plddt.npy", ".meta.json"):
        fp = f"{outdir}/{tag}{ext}"
        if not os.path.exists(fp) or os.path.getsize(fp) == 0:
            return False
    return True


def _write_outputs(outdir, tag, pdb, pae, pl, pairs, labels, vid, rc, lk, order, t0):
    open(f"{outdir}/{tag}.pdb", "w", encoding="utf-8").write(pdb)
    np.save(f"{outdir}/{tag}.pae.npy", pae)
    np.save(f"{outdir}/{tag}.plddt.npy", pl)
    meta = {"labels": labels, "chain_ids": [c for c, _ in pairs],
            "seqs": [s for _, s in pairs], "lens": [len(s) for _, s in pairs],
            "binder_chain": pairs[0][0], "variant": vid, "recycles": rc,
            "linker": lk, "chain_order": order, "mean_plddt": float(pl.mean()),
            "pae_shape": list(pae.shape), "seconds": round(time.time() - t0, 1)}
    json.dump(meta, open(f"{outdir}/{tag}.meta.json", "w", encoding="utf-8"))
    tot = sum(len(s) for _, s in pairs)
    ok = pae.shape[0] == len(labels) == tot
    print(f"[{'OK ' if ok else '⚠ 長度不符'}] {tag}  pLDDT {meta['mean_plddt']:.1f}  "
          f"{meta['seconds']:.1f}s  PAE {pae.shape}  鏈長和 {tot}", flush=True)
    return meta


def cmd_fold(a):
    import torch
    _seed_warning(a)
    a.device = pick_device(a.device)
    print(f"device = {a.device}   recycles={a.recycles}  linkers={a.linkers}  "
          f"orders={a.chain_orders}", flush=True)
    os.makedirs(a.out, exist_ok=True)
    jobs = read_spec(a.spec)
    idx, first = [], True
    pending_tasks = []

    for name, pairs in jobs:
        for lk in a.linkers:
            for order in a.chain_orders:
                feed = pairs if order == "fwd" else list(reversed(pairs))
                full, pos, lab = build(feed, lk)
                keep = np.array([i for i, l in enumerate(lab) if l is not None])
                labels = [lab[i] for i in keep]
                for rc in a.recycles:
                    vid = f"r{rc}L{lk}" + ("rev" if order == "rev" else "")
                    tag = f"{name}__{vid}"
                    if not getattr(a, "overwrite", False) and _is_tag_completed(a.out, tag):
                        try:
                            meta = json.load(open(f"{a.out}/{tag}.meta.json", encoding="utf-8"))
                            print(f"[SKIP 斷點續跑] {tag} 已存在 (pLDDT {meta.get('mean_plddt', 0):.1f})，直接沿用", flush=True)
                            idx.append({"name": name, "tag": tag, **meta})
                            continue
                        except Exception:
                            pass
                    pending_tasks.append({
                        "name": name, "pairs": pairs, "lk": lk, "order": order,
                        "rc": rc, "vid": vid, "tag": tag, "full": full, "pos": pos,
                        "lab": lab, "keep": keep, "labels": labels,
                    })

    if not pending_tasks:
        json.dump(idx, open(f"{a.out}/index.json", "w", encoding="utf-8"), indent=1)
        return

    tok, model = load_hf(a.device, a.chunk, a.fp16)

    # 將待算任務按 (rc, 序列總長) 排序聚合，讓同長度的不同抗體自動組成 Batch 並行推理
    if len(pending_tasks) > 1 and not a.check_determinism:
        pending_tasks.sort(key=lambda x: (x["rc"], len(x["full"])))

    # 依 (rc, 序列總長) 分組進行同長度動態 Batch 並行推理，壓榨 A100 顯存與 Tensor Cores
    i = 0
    while i < len(pending_tasks):
        t_head = pending_tasks[i]
        bmax = 1 if (a.check_determinism and first) else _auto_batch_size(a.device, len(t_head["full"]))
        batch = [t_head]
        j = i + 1
        while len(batch) < bmax and j < len(pending_tasks):
            cand = pending_tasks[j]
            if cand["rc"] == t_head["rc"] and len(cand["full"]) == len(t_head["full"]):
                batch.append(cand)
                j += 1
            else:
                break
        i = j
        t0 = time.time()
        full_list = [x["full"] for x in batch]
        pos_list = [x["pos"] for x in batch]
        rc = t_head["rc"]

        if len(batch) > 1:
            print(f"[A100 Batch={len(batch)}] 並行推理 {[x['tag'] for x in batch]} (len={len(full_list[0])}, rc={rc})...", flush=True)

        try:
            batch_outs = infer_hf_batch(model, tok, full_list, pos_list, a.device, rc)
        except Exception as e:
            # 若 Batch 推理 OOM，自動退回逐條推理
            if a.device == "cuda":
                torch.cuda.empty_cache()
            batch_outs = []
            for x in batch:
                try:
                    batch_outs.append(infer_hf(model, tok, x["full"], x["pos"], a.device, rc))
                except Exception as e2:
                    if "out of memory" not in str(e2).lower():
                        print(f"[FAIL] {x['tag']}: {type(e2).__name__}: {e2}", flush=True)
                        batch_outs.append(None)
                        continue
                    if a.device == "cuda":
                        torch.cuda.empty_cache()
                    small = max(8, a.chunk // 4)
                    print(f"[OOM] {x['tag']}: 改用 chunk={small} 重試", flush=True)
                    model.trunk.set_chunk_size(small)
                    try:
                        batch_outs.append(infer_hf(model, tok, x["full"], x["pos"], a.device, rc))
                    except Exception as e3:
                        print(f"[FAIL] {x['tag']}: {type(e3).__name__}: {e3}", flush=True)
                        batch_outs.append(None)
                    finally:
                        model.trunk.set_chunk_size(a.chunk or None)

        if a.check_determinism and first and batch_outs and batch_outs[0] is not None:
            pae2 = infer_hf(model, tok, batch[0]["full"], batch[0]["pos"], a.device, rc)[0]
            print(f"[determinism] 同設定跑兩次 PAE 最大差 "
                  f"{np.abs(batch_outs[0][0] - pae2).max():.3e}  "
                  f"（≈0 即證明 seed 無效，需用變體測穩定度）", flush=True)
            first = False

        for x, out_item in zip(batch, batch_outs):
            if out_item is None:
                continue
            pae, pl, raw = out_item
            keep = x["keep"]
            pae, pl = pae[np.ix_(keep, keep)], pl[keep]
            try:
                pdb = clean_pdb(raw, x["lab"])
            except ValueError as e:
                print(f"[FAIL] {x['tag']}: PDB 清理失敗：{e}", flush=True)
                continue
            meta = _write_outputs(a.out, x["tag"], pdb, pae, pl, x["pairs"], x["labels"],
                                  x["vid"], x["rc"], x["lk"], x["order"], t0)
            idx.append({"name": x["name"], "tag": x["tag"], **meta})
            json.dump(idx, open(f"{a.out}/index.json", "w", encoding="utf-8"), indent=1)

    json.dump(idx, open(f"{a.out}/index.json", "w", encoding="utf-8"), indent=1)


# ═══════════════════════════════════════════════════════════════════
#  6. 摺疊：fair-esm 後端（ModelScope / 24GB 卡）
# ═══════════════════════════════════════════════════════════════════
# 官方配置：1 seed、1 sample、無模板、僅靶點序列用 MSA、設計不用 MSA。
# ESMFold 本身不吃 MSA（用 ESM-2 表徵），故「MSA」一項無對應操作 ——
# 此為已去函主辦方確認的項目之一，若官方另有前處理需在此補上。

def cmd_fold_esm(a):
    import torch, esm
    _seed_warning(a)
    os.makedirs(a.out, exist_ok=True)
    jobs = read_spec(a.spec)
    model = esm.pretrained.esmfold_v1().eval()
    if a.device == "cuda":
        model = model.cuda()
    if a.chunk:
        model.set_chunk_size(a.chunk)

    idx = []
    for name, pairs in jobs:
        for lk in a.linkers:
            for order in a.chain_orders:
                feed = pairs if order == "fwd" else list(reversed(pairs))
                feed_labels = [c for c, _ in feed]
                for rc in a.recycles:
                    vid = f"r{rc if rc is not None else 'def'}L{lk}" + \
                          ("rev" if order == "rev" else "")
                    tag, t0 = f"{name}__{vid}", time.time()
                    try:
                        with torch.no_grad():
                            out = model.infer(sequences=[":".join(s for _, s in feed)],
                                              num_recycles=rc,
                                              chain_linker="G" * lk)
                        raw = model.output_to_pdb(out)[0]
                        pae = out["predicted_aligned_error"][0].cpu().numpy()
                        pl = plddt_per_res(out["plddt"][0].cpu().numpy())
                        if "chain_index" in out:
                            ci = np.asarray(out["chain_index"][0].cpu().numpy())
                        else:                       # 退路：依序列長度推回
                            ci = []
                            for k, (_, s) in enumerate(feed):
                                ci += [k] * len(s)
                                if k < len(feed) - 1:
                                    ci += [-1] * lk
                            ci = np.array(ci)
                    except Exception as e:
                        print(f"[FAIL] {tag}: {type(e).__name__}: {e}"); continue

                    lab = [None if int(c) < 0 else feed_labels[int(c)] for c in ci]
                    keep = np.array([i for i, l in enumerate(lab) if l is not None])
                    if keep.size == 0 or pae.shape[0] != len(lab):
                        print(f"[FAIL] {tag}: chain_index 長度 {len(lab)} 與 "
                              f"PAE {pae.shape} 不符"); continue
                    labels = [lab[i] for i in keep]
                    pae, pl = pae[np.ix_(keep, keep)], pl[keep]
                    # fair-esm 的 output_to_pdb 行為若與假設不符，保留原始檔並明確警告，
                    # 不默默輸出錯編號。
                    try:
                        pdb = clean_pdb(raw, lab)
                    except ValueError as e:
                        pdb = raw
                        print(f"⚠ {tag}: PDB 未清理（{e}），疊合前請先確認鏈與編號")
                    meta = _write_outputs(a.out, tag, pdb, pae, pl, pairs, labels,
                                          vid, rc, lk, order, t0)
                    idx.append({"name": name, "tag": tag, **meta})
    json.dump(idx, open(f"{a.out}/index.json", "w", encoding="utf-8"), indent=1)


# ═══════════════════════════════════════════════════════════════════
#  7. 匯總計分
# ═══════════════════════════════════════════════════════════════════

def _mstd(vals):
    v = np.asarray([x for x in vals if x is not None], dtype=float)
    return "—" if v.size == 0 else f"{v.mean():.3f} ± {v.std():.3f}"


def cmd_score(a):
    main_cut, alt_cuts = a.cutoffs[0], a.cutoffs[1:]
    rows = []
    for f in sorted(glob.glob(f"{a.dir}/*.pae.npy")):
        tag = os.path.basename(f)[:-8]
        mp = f"{a.dir}/{tag}.meta.json"
        if not os.path.exists(mp):
            print(f"跳過 {tag}: 缺 meta.json"); continue
        meta = json.load(open(mp, encoding="utf-8"))
        labels = meta["labels"]
        binder = a.binder or meta.get("binder_chain") or sorted(set(labels))[0]
        pae = np.load(f)
        if pae.shape[0] != len(labels):
            print(f"跳過 {tag}: PAE {pae.shape} 與標籤數 {len(labels)} 不符"); continue
        r = {"name": tag.split("__")[0], "tag": tag, "variant": meta.get("variant", ""),
             "recycles": meta.get("recycles"), "linker": meta.get("linker"),
             "chain_order": meta.get("chain_order"), "binder": binder,
             "mean_plddt": meta.get("mean_plddt"), "n_chains": len(set(labels))}
        if len(set(labels)) < 2:
            r["pair"] = "—"               # 單鏈靶點沒有鏈間界面可計分
        else:
            for c in a.cutoffs:
                s = summarize(ipsae(pae, labels, c, binder=binder), binder)
                pre = f"ipsae{c:g}"
                r[f"{pre}_max"], r[f"{pre}_min"] = s["pair_max"], s["pair_min"]
                r[f"{pre}_n"] = s["n_pass"]
                r[f"{pre}_merged_max"], r[f"{pre}_merged_n"] = s["merged_max"], s["merged_n"]
                if c == main_cut:
                    r["pair"] = s["pair"]
                    r["all_pairs"] = " ".join(f"{k}={v:.3f}"
                                              for k, v in s["all_pairs"].items())
        pdb = f"{a.dir}/{tag}.pdb"
        r["iface"] = None
        if os.path.exists(pdb) and len(set(labels)) >= 2:
            try:
                r["iface"] = interface(pdb, binder, contact_cutoff=a.contact_cutoff)
            except Exception as e:
                print(f"⚠ {tag}: 界面指標失敗 {type(e).__name__}: {e}")
        seqs, cids = meta.get("seqs"), meta.get("chain_ids")
        r["triage"] = (triage(seqs[cids.index(binder)], a.max_len)
                       if seqs and cids and binder in cids else None)
        rows.append(r)

    if not rows:
        print("沒有可計分的輸出。"); return
    groups = {}
    for r in rows:
        groups.setdefault(r["name"], []).append(r)
    order = sorted(groups, key=lambda n: -np.mean(
        [g.get(f"ipsae{main_cut:g}_max") or 0.0 for g in groups[n]]))

    print(f"表一 計分 —— 主 cutoff {main_cut:g} Å"
          + (f"，對照 {', '.join(f'{c:g}' for c in alt_cuts)} Å" if alt_cuts else "")
          + "（不同 cutoff 的數值不可互比）\n")
    hdr = f"{'名稱':<30}{'n':>3}{'鏈對':>7}{'ipSAE_max':>20}{'ipSAE_min':>20}"
    for c in alt_cuts:
        hdr += f"{f'max@{c:g}':>20}"
    print(hdr + f"{'合併A-ALL':>20}")
    pm = f"ipsae{main_cut:g}"
    for n in order:
        g = groups[n]
        line = (f"{n:<30}{len(g):>3}{(g[0].get('pair') or '—'):>7}"
                f"{_mstd([x.get(pm+'_max') for x in g]):>20}"
                f"{_mstd([x.get(pm+'_min') for x in g]):>20}")
        for c in alt_cuts:
            line += f"{_mstd([x.get(f'ipsae{c:g}_max') for x in g]):>20}"
        print(line + f"{_mstd([x.get(pm+'_merged_max') for x in g]):>20}")

    multi = [n for n in order if groups[n][0]["n_chains"] > 2]
    if multi:
        print("\n多聚體靶點的各鏈對明細（主 cutoff，第一個變體）：")
        for n in multi:
            print(f"  {n:<30}{groups[n][0].get('all_pairs', '')}")

    print(f"\n表二 界面幾何與序列 triage（接觸定義：CB–CB < {a.contact_cutoff:g} Å；"
          f"長度上限 {a.max_len} aa）\n")
    print(f"{'名稱':<30}{'接觸':>6}{'min重原子':>11}{'pDockQ':>9}{'pLDDT':>8}"
          f"{'binder長':>9}{'C數':>5}{'淨電荷':>8}{'GRAVY':>8}  旗標")
    for n in order:
        g = groups[n]
        ifs = [x["iface"] for x in g if x["iface"]]
        t = next((x["triage"] for x in g if x["triage"]), None)
        con = f"{np.mean([i['n_contacts_cb8'] for i in ifs]):.0f}" if ifs else "—"
        mh = f"{np.mean([i['min_heavy_dist'] for i in ifs]):.2f}" if ifs else "—"
        pdq = f"{np.mean([i['pdockq'] for i in ifs]):.3f}" if ifs else "—"
        pls = [x["mean_plddt"] for x in g if x["mean_plddt"] is not None]
        pl = f"{np.mean(pls):.1f}" if pls else "—"
        if t:
            blen, ncys = str(t["len"]), str(t["n_cys"])
            chg, gv = f"{t['net_charge']:+.1f}", f"{t['gravy']:+.2f}"
            fl = "、".join(t["flags"]) if t["flags"] else "OK"
        else:
            blen = ncys = chg = gv = fl = "—"
        print(f"{n:<30}{con:>6}{mh:>11}{pdq:>9}{pl:>8}"
              f"{blen:>9}{ncys:>5}{chg:>8}{gv:>8}  {fl}")

    print("\n判讀提示：")
    print("  · 接觸數 0 或 min 重原子 > 5 Å → 兩條鏈根本沒對上，ipSAE 低不是 cutoff 的問題。")
    print("  · pDockQ > 0.23 約為可接受界面、> 0.5 多為正確界面（Bryant 2022）。")
    print("  · 合併 A-ALL 明顯高於 pair_max → 跨亞基結合，按鏈對計分會低估。")
    print("  · 旗標非 OK 的序列應在投稿前先修（長度／未配對 C／疏水聚集）。")

    if a.csv:
        flat = []
        for r in rows:
            d = {k: v for k, v in r.items() if k not in ("iface", "triage")}
            d.update({f"iface_{k}": v for k, v in (r["iface"] or {}).items()})
            tr = dict(r["triage"] or {})
            if "flags" in tr:
                tr["flags"] = "|".join(tr["flags"])
            d.update({f"triage_{k}": v for k, v in tr.items()})
            flat.append(d)
        cols = sorted({k for d in flat for k in d})
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(flat)
        print(f"\nCSV 已寫出：{a.csv}（{len(flat)} 列 × {len(cols)} 欄）")


# ═══════════════════════════════════════════════════════════════════
#  8. 與晶體結構疊合
# ═══════════════════════════════════════════════════════════════════
SITES = {"位點 I primary": [73, 75, 77, 88, 89, 90, 91, 92, 97, 135, 136, 137, 138, 139, 140],
         "位點 II": [20, 21, 22, 23, 24, 25, 65, 66, 67, 70, 105, 107, 111, 137, 138, 139, 140, 141, 142],
         "位點 I 熱點": [90, 92, 97, 135, 136, 137],
         "位點 II 熱點": [23, 24, 65, 138, 140, 141]}


def fetch_pdb(pdb_id, dest=None):
    pdb_id = pdb_id.upper()
    dest = dest or f"{pdb_id}.pdb"
    if os.path.exists(dest):
        print(f"沿用既有檔案 {dest}")
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"下載 {url} → {dest}")
        urllib.request.urlretrieve(url, dest)
    return dest


def ca_chains(path):
    """{chain: {resSeq: (aa, CA 座標)}}"""
    out = {}
    for ch, rs in read_pdb(path).items():
        d = {r["resseq"]: (r["aa"], r["ca"]) for r in rs if r["ca"] is not None}
        if d:
            out[ch] = d
    return out


def best_offset(pred, ref):
    """(offset, 相符數)：使 pred[k] 與 ref[k+offset] 殘基種類相符者最多。
    晶體檔缺殘基或用 UniProt 編號時都不會整段錯位。"""
    votes = {}
    for k, (aa, _) in pred.items():
        for j, (bb, _) in ref.items():
            if aa == bb:
                votes[j - k] = votes.get(j - k, 0) + 1
    if not votes:
        return 0, 0
    off = max(votes, key=lambda o: votes[o])
    return off, votes[off]


def matched_keys(pred, ref, offset):
    return [k for k in sorted(pred)
            if (k + offset) in ref and pred[k][0] == ref[k + offset][0]]


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, W = np.linalg.svd(Pc.T @ Qc)
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ W))]) @ W
    return Pc @ R, Qc, R


def cmd_compare(a):
    ref_path = a.ref or (fetch_pdb(a.fetch) if a.fetch else None)
    if not ref_path:
        sys.exit("請給 ref 檔路徑，或用 --fetch 5M2J 自動下載")
    pc, rc = ca_chains(a.pred), ca_chains(ref_path)
    if not pc or not rc:
        sys.exit("其中一個檔案沒讀到 CA 原子")

    combos = []
    for pk in ([a.pred_chain] if a.pred_chain != "auto" else sorted(pc)):
        for rk in ([a.ref_chain] if a.ref_chain != "auto" else sorted(rc)):
            if pk not in pc or rk not in rc:
                continue
            pred = {k - a.pred_offset: v for k, v in pc[pk].items()}
            auto = a.auto_offset and a.ref_offset is None
            off = best_offset(pred, rc[rk])[0] if auto else -(a.ref_offset or 0)
            combos.append((len(matched_keys(pred, rc[rk], off)), pk, rk, off, pred))
    if not combos:
        sys.exit("指定的鏈不存在")
    nmatch, pk, rk, off, pred = max(combos, key=lambda c: c[0])
    ref = rc[rk]
    keys = matched_keys(pred, ref, off)

    print(f"預測 {os.path.basename(a.pred)} 鏈 {pk}（{len(pc[pk])} 殘基，"
          f"共 {len(pc)} 條鏈：{''.join(sorted(pc))}）")
    print(f"參考 {os.path.basename(ref_path)} 鏈 {rk}（{len(ref)} 殘基，"
          f"共 {len(rc)} 條鏈：{''.join(sorted(rc))}）")
    auto = a.auto_offset and a.ref_offset is None
    print(f"編號偏移：ref = pred {off:+d}   殘基種類相符 {len(keys)} 個"
          f"{'（自動對位）' if auto else '（手動指定）'}")
    if len(keys) < 20:
        sys.exit(f"相符殘基僅 {len(keys)} 個，鏈或編號可能不匹配"
                 f"（試 --ref-chain / --ref-offset）")

    P = np.array([pred[k][1] for k in keys])
    Q = np.array([ref[k + off][1] for k in keys])
    Ac, Bc = P.mean(0), Q.mean(0)
    Pr, Qc_, R = kabsch(P, Q)
    print(f"整體 CA RMSD（全鏈疊合）: "
          f"{float(np.sqrt(((Pr - Qc_) ** 2).sum() / len(keys))):.2f} Å\n")

    print(f"  {'區域':<18}{'殘基':>5}{'RMSD':>9}")
    for nm, res in SITES.items():
        ks = [k for k in res if k in keys]
        if not ks:
            print(f"  {nm:<18}{0:>5}{'—':>9}"); continue
        Pp = np.array([pred[k][1] for k in ks]) - Ac
        Qq = np.array([ref[k + off][1] for k in ks]) - Bc
        r = float(np.sqrt((((Pp @ R) - Qq) ** 2).sum() / len(ks)))
        print(f"  {nm:<18}{len(ks):>5}{f'{r:.2f} Å':>9}")

    dev = sorted(((float(np.linalg.norm((pred[k][1] - Ac) @ R - (ref[k + off][1] - Bc))), k)
                  for k in keys), reverse=True)
    print("\n偏移最大的 12 個殘基（預測檔編號 / 偏移 Å）：")
    print("  " + ", ".join(f"{k}:{d:.1f}" for d, k in dev[:12]))
    print("\n判讀：表位區 RMSD < 1.0 Å → 可用晶體結構設計；"
          "> 1.5 Å → 應改用預測結構作為設計靶點。")


# ═══════════════════════════════════════════════════════════════════
#  9. 自測（不需 GPU）
# ═══════════════════════════════════════════════════════════════════

def _ref_directional(pae, idx_from, idx_to, cutoff):
    """逐列迴圈的參考實作，僅用於驗證向量化版本一致。"""
    best = 0.0
    for i in idx_from:
        row = np.asarray(pae)[i, idx_to].astype(float)
        keep = row < cutoff
        n = int(keep.sum())
        if n:
            best = max(best, float(np.mean(1.0 / (1.0 + (row[keep] / d0_from_n(n)) ** 2))))
    return best


def _fake_hf_pdb(full, pos, plddt=80.0):
    """模擬 HF output_to_pdb：單一鏈 A、resSeq 直接用 position_ids+1（含 +512 跳號）。"""
    three = {v: k for k, v in AA3.items() if k != "MSE"}
    lines, ser = [], 0
    for i, aa in enumerate(full):
        for an, dx, el in (("N", 0.0, "N"), ("CA", 1.0, "C"), ("CB", 2.0, "C")):
            ser += 1
            lines.append(f"ATOM  {ser:5d}  {an:<3s}{' '}{three[aa]} A{pos[i]+1:4d}    "
                         f"{i*3.8+dx:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{plddt:6.2f}"
                         f"           {el}")
    lines += [f"TER   {ser+1:5d}      {three[full[-1]]} A{pos[-1]+1:4d}", "END"]
    return "\n".join(lines) + "\n"


def cmd_selftest(a):
    ok = True

    def chk(label, got, exp, tol=0.02):
        nonlocal ok
        good = abs(got - exp) <= tol
        ok &= good
        print(f"  [{'OK  ' if good else 'FAIL'}] {label}: 得 {got:.3f}，期望 ≈{exp:.3f}")

    def flag(label, good):
        nonlocal ok
        ok &= bool(good)
        print(f"  [{'OK  ' if good else 'FAIL'}] {label}")

    print("性質 1 — 小界面被懲罰（文獻基準：PAE 2 Å 時 n=50 得 ≈0.56、n=600 得 ≈0.95）")
    for n, exp in ((50, 0.56), (600, 0.95)):
        pae = np.full((n + 1, n + 1), 2.0)
        chk(f"n={n}", ipsae(pae, ["A"] + ["B"] * n, 10.0)["A-B"]["asym_a_to_b"], exp)

    print("\n性質 2 — 無序／不接觸區不拉低分數（ipTM 會，ipSAE 不會）")
    n = 50; L = n + 1
    pae = np.full((L, L), 2.0)
    chain = ["A"] + ["B"] * n
    base = ipsae(pae, chain, 10.0)["A-B"]["asym_a_to_b"]
    pae2 = np.full((L + 120, L + 120), 30.0); pae2[:L, :L] = pae
    aug = ipsae(pae2, chain + ["B"] * 120, 10.0)["A-B"]["asym_a_to_b"]
    chk("加入 120 個無序殘基後分數不變", aug, base, tol=1e-6)

    print("\n性質 3 — 非對稱：A→B 與 B→A 可不同")
    pae = np.full((41, 41), 25.0)
    pae[0, 1:41] = 1.5; pae[1:41, 0] = 9.0
    r = ipsae(pae, ["A"] + ["B"] * 40, 10.0)["A-B"]
    print(f"  A→B {r['asym_a_to_b']:.3f}   B→A {r['asym_b_to_a']:.3f}   "
          f"max {r['max']:.3f}  min {r['min']:.3f}")
    flag(f"兩方向差異 {abs(r['asym_a_to_b']-r['asym_b_to_a']):.3f} > 0.1",
         abs(r["asym_a_to_b"] - r["asym_b_to_a"]) > 0.1)

    print("\n性質 4 — cutoff 改變會改變數值（故不同 cutoff 不可互比）")
    pae = np.full((101, 101), 2.0); pae[0, 51:] = 12.0
    chain = ["A"] + ["B"] * 100
    r10 = ipsae(pae, chain, 10.0)["A-B"]["asym_a_to_b"]
    r15 = ipsae(pae, chain, 15.0)["A-B"]["asym_a_to_b"]
    print(f"  cutoff 10 Å → {r10:.3f}   cutoff 15 Å → {r15:.3f}   差 {abs(r10-r15):.3f}")
    flag("兩 cutoff 給出不同數值", abs(r10 - r15) > 0.01)

    print("\n性質 5 — 向量化與逐列迴圈參考實作完全一致")
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(20):
        L = int(rng.integers(30, 90))
        pae = rng.uniform(0.3, 30.0, size=(L, L))
        cut = int(L * rng.uniform(0.2, 0.8))
        ia, ib = np.arange(cut), np.arange(cut, L)
        worst = max(worst, abs(ipsae_directional(pae, ia, ib, 10.0)[0]
                               - _ref_directional(pae, ia, ib, 10.0)))
    flag(f"20 組隨機矩陣最大誤差 {worst:.2e}", worst < 1e-12)

    print("\n性質 6 — 跨亞基結合：pairwise 會低估，合併靶點讀數較高")
    L = 10 + 40 * 3
    pae = np.full((L, L), 30.0)
    chain = ["A"] * 10 + ["B"] * 40 + ["C"] * 40 + ["D"] * 40
    pae[0, 10:30] = 2.0; pae[10:30, 0] = 2.0
    pae[0, 50:70] = 2.0; pae[50:70, 0] = 2.0
    s = summarize(ipsae(pae, chain, 10.0), "A")
    print("  各鏈對 max: " + "  ".join(f"{k} {v:.3f}" for k, v in s["all_pairs"].items()))
    print(f"  pair_max {s['pair_max']:.3f}（n={s['n_pass']}）   "
          f"合併 A-BCD {s['merged_max']:.3f}（n={s['merged_n']}）")
    flag("合併讀數高於最佳鏈對，且 pair_max 不是只看 A-B",
         s["merged_max"] > s["pair_max"] + 0.05 and s["pair"] in ("A-B", "A-C"))

    print("\n性質 7 — 單鏈輸入回傳空字典（無鏈間界面可計分）")
    flag("得 {}", ipsae(np.full((20, 20), 2.0), ["A"] * 20, 10.0) == {})

    print("\n性質 8 — PDB 清理：剔除 linker、分鏈、各鏈殘基編號從 1 起算")
    pairs = [("A", "MKVWY"), ("B", "DEFIL"), ("C", "PSTQN")]
    full, pos, lab = build(pairs, linker_len=4)
    raw = _fake_hf_pdb(full, pos)
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "c.pdb")
        open(p, "w", encoding="utf-8").write(clean_pdb(raw, lab))
        ch = read_pdb(p)
        got = {c: "".join(r["aa"] for r in rs) for c, rs in sorted(ch.items())}
        nums = {c: (rs[0]["resseq"], rs[-1]["resseq"]) for c, rs in sorted(ch.items())}
        print(f"  原始 resSeq 含 +{JUMP} 跳號：{pos[0]+1}…{pos[-1]+1}（{len(full)} 殘基）")
        print(f"  清理後 {got}  編號 {nums}")
        flag("三條鏈、序列正確、各鏈從 1 起算",
             got == {"A": "MKVWY", "B": "DEFIL", "C": "PSTQN"}
             and all(v == (1, 5) for v in nums.values()))
        iface_ok = interface(p, "A") is not None
        flag("介面指標可在清理後的檔案上執行", iface_ok)
    try:
        clean_pdb(raw, lab[:-1]); flag("殘基數不符時應該擋下來", False)
    except ValueError:
        flag("殘基數不符時擋下來，不默默輸出錯編號", True)

    print("\n性質 9 — triage 攔截長度超標與未配對 Cysteine")
    flag("長度 260 被擋", "長度260>250" in "".join(triage("A" * 260)["flags"]))
    flag("奇數 C 被擋", triage("MKVWYC")["cys_odd"])
    flag("正常奈米抗體序列過關", triage("QVQLQESGGGLVQPGGSLRLSCAASGRTFSDHSGYTYT"
                                        "IGWFRQAPGKEREFVARIYWSSGNTYYADSVKGRFAIS"
                                        "RDIAKNTVDLTMNNLEPEDTAVYYCAARDGIPTSRSVE"
                                        "SYNYWGQGTQVTVSS")["pass"])

    print(f"\n自測整體：{'全部通過' if ok else '有失敗項'}")
    return ok


# ═══════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════

def _add_fold_args(sp):
    sp.add_argument("--spec", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--linkers", type=int, nargs="+", default=[25],
                    help="G-linker 長度變體，例如 --linkers 20 25 30")
    sp.add_argument("--chain-orders", nargs="+", default=["fwd"], choices=["fwd", "rev"],
                    help="鏈輸入順序；鏈 ID 仍依 spec 固定")
    sp.add_argument("--seeds", type=int, nargs="+", default=None,
                    help="保留給舊指令；ESMFold 推理是確定性的，多值會摺疊成一次")
    sp.add_argument("--chunk", type=int, default=64,
                    help="trunk 分塊大小，省顯存但變慢；0 = 不分塊（40 GB 以上建議）")
    sp.add_argument("--overwrite", action="store_true",
                    help="強制重算已存在的輸出檔（預設自動跳過已完成的序列變體以支援斷點續跑）")


def main():
    p = argparse.ArgumentParser(
        description="TNF-α ipSAE 測試管線（單檔版）",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("selftest", help="不需 GPU 的自測").set_defaults(fn=cmd_selftest)

    sp = sub.add_parser("triage", help="序列預檢")
    sp.add_argument("--spec", nargs="*", default=None)
    sp.add_argument("--seq", default=None)
    sp.add_argument("--max-len", type=int, default=250)
    sp.set_defaults(fn=cmd_triage)

    sp = sub.add_parser("fold", help="ESMFold 推理（HuggingFace 後端）")
    _add_fold_args(sp)
    sp.add_argument("--recycles", type=int, nargs="+", default=[4],
                    help="回收次數變體；預設 4 = HF 模型預設值")
    sp.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    sp.add_argument("--fp16", dest="fp16", action="store_true", default=True)
    sp.add_argument("--no-fp16", dest="fp16", action="store_false")
    sp.add_argument("--check-determinism", action="store_true",
                    help="第一個任務跑兩次，印出兩次 PAE 的最大差異")
    sp.set_defaults(fn=cmd_fold)

    sp = sub.add_parser("fold-esm", help="ESMFold 推理（fair-esm 後端）")
    _add_fold_args(sp)
    sp.add_argument("--recycles", type=int, nargs="+", default=[None],
                    help="回收次數變體；預設 None = 模型預設值")
    sp.add_argument("--device", default="cuda")
    sp.set_defaults(fn=cmd_fold_esm)

    sp = sub.add_parser("score", help="匯總計分")
    sp.add_argument("--dir", required=True)
    sp.add_argument("--cutoffs", type=float, nargs="+", default=[10.0, 15.0],
                    help="第一個為主 cutoff，其餘為對照")
    sp.add_argument("--binder", default=None, help="強制指定 binder 鏈 ID")
    sp.add_argument("--max-len", type=int, default=250)
    sp.add_argument("--contact-cutoff", type=float, default=8.0)
    sp.add_argument("--csv", default=None)
    sp.set_defaults(fn=cmd_score)

    sp = sub.add_parser("compare", help="與晶體結構疊合")
    sp.add_argument("pred")
    sp.add_argument("ref", nargs="?", default=None)
    sp.add_argument("--fetch", default=None, help="RCSB ID，例如 5M2J / 1TNF / 2AZ5")
    sp.add_argument("--pred-chain", default="auto")
    sp.add_argument("--ref-chain", default="auto")
    sp.add_argument("--pred-offset", type=int, default=0)
    sp.add_argument("--ref-offset", type=int, default=None,
                    help="手動指定編號偏移（給定時不做自動對位）")
    sp.add_argument("--no-auto-offset", dest="auto_offset",
                    action="store_false", default=True)
    sp.set_defaults(fn=cmd_compare)

    a = p.parse_args()
    r = a.fn(a)
    if a.cmd == "selftest":
        raise SystemExit(0 if r else 1)


def run_cli(argv):
    """供 worker.py 在同一進程內直接調用，讓 ESMFold 常駐顯存免重複載入。"""
    old_argv = sys.argv
    try:
        sys.argv = ["tnf_pipeline.py"] + list(argv)
        try:
            main()
            return True
        except SystemExit as e:
            return e.code in (0, None)
    finally:
        sys.argv = old_argv


if __name__ == "__main__":
    main()
