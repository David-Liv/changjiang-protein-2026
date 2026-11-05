"""純標準庫 PDB 解析（容器無 biopython／PyPI 不可達）。固定欄位解析，不做任何猜測。"""
import math, re
from collections import OrderedDict, defaultdict

class Atom:
    __slots__=("rec","serial","name","altloc","resname","chain","resseq","icode",
               "x","y","z","occ","b","element")
    def __init__(s,l):
        s.rec=l[0:6].strip(); s.serial=int(l[6:11]); s.name=l[12:16].strip()
        s.altloc=l[16]; s.resname=l[17:20].strip(); s.chain=l[21]
        s.resseq=int(l[22:26]); s.icode=l[26]
        s.x=float(l[30:38]); s.y=float(l[38:46]); s.z=float(l[46:54])
        s.occ=float(l[54:60] or 1.0); s.b=float(l[60:66] or 0.0)
        s.element=l[76:78].strip() or s.name[0]
    @property
    def key(s): return (s.chain, s.resseq, s.icode)
    @property
    def is_h(s): return s.element == "H"

AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split())
THREE2ONE = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split(),
                     list("ARNDCQEGHILKMFPSTWYVM")))

class Structure:
    def __init__(self, path, keep_altloc_first=True):
        self.path=path; self.atoms=[]; self.remark465=[]; self.resolution=None
        self.compnd=[]; self.het_names={}; self.links=[]; self.seqadv=[]; self.dbref=[]
        seen=set()
        with open(path) as fh:
            for l in fh:
                t=l[:6]
                if t in ("ATOM  ","HETATM"):
                    if len(l)<54: continue
                    a=Atom(l)
                    if a.is_h: continue
                    if keep_altloc_first:
                        k=(a.chain,a.resseq,a.icode,a.name)
                        if a.altloc not in (" ","A") : continue
                        if k in seen: continue
                        seen.add(k)
                    self.atoms.append(a)
                elif l.startswith("REMARK 465"):  self.remark465.append(l.rstrip())
                elif l.startswith("REMARK   2 RESOLUTION"):
                    m=re.search(r"RESOLUTION\.\s*([0-9.]+)\s*ANGSTROM", l)
                    if m: self.resolution=float(m.group(1))
                elif l.startswith("COMPND"): self.compnd.append(l[10:].rstrip())
                elif l.startswith("HETNAM"): self.het_names[l[11:14].strip()]=l[15:].strip()
                elif l.startswith("LINK"):   self.links.append(l.rstrip())
                elif l.startswith("SEQADV"): self.seqadv.append(l.rstrip())
                elif l.startswith("DBREF"):  self.dbref.append(l.rstrip())
        self.by_chain=defaultdict(list)
        for a in self.atoms: self.by_chain[a.chain].append(a)

    def residues(self, chain=None, protein_only=False):
        """OrderedDict[(chain,resseq,icode)] -> (resname, [atoms])"""
        out=OrderedDict()
        for a in self.atoms:
            if chain and a.chain!=chain: continue
            if protein_only and a.resname not in AA3: continue
            out.setdefault(a.key,(a.resname,[]))[1].append(a)
        return out

    def missing_residues(self):
        """從 REMARK 465 取 (resname, chain, resseq)。只讀表格列。"""
        out=[]; started=False
        for l in self.remark465:
            if "M RES C SSSEQI" in l: started=True; continue
            if not started: continue
            f=l[15:].split()
            if len(f)>=3 and f[0].isalpha() and len(f[0])==3:
                try: out.append((f[0], f[1], int(f[2])))
                except ValueError: pass
            elif len(f)>=4 and f[1].isalpha():
                try: out.append((f[1], f[2], int(f[3])))
                except ValueError: pass
        return out

def dist2(a,b): return (a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2

class Grid:
    """均勻網格鄰居搜尋，避免 O(N²)。"""
    def __init__(self, atoms, cell):
        self.cell=cell; self.g=defaultdict(list)
        for a in atoms: self.g[(int(a.x//cell),int(a.y//cell),int(a.z//cell))].append(a)
    def near(self, a, cutoff):
        c=self.cell; r=int(math.ceil(cutoff/c))
        i0,j0,k0=int(a.x//c),int(a.y//c),int(a.z//c); cut2=cutoff*cutoff; out=[]
        for i in range(i0-r,i0+r+1):
            for j in range(j0-r,j0+r+1):
                for k in range(k0-r,k0+r+1):
                    for b in self.g.get((i,j,k),()):
                        if dist2(a,b)<=cut2: out.append(b)
        return out
