import sys, types
 
# Windows'un Uygulama Denetimi politikasi, Bio.PDB'nin otomatik yuklemeye
# calistigi Bio.Align alt modulunu (biz hic kullanmiyoruz) engelliyor.
# Bunu atlatmak icin, gercek modul yuklenmeden once sahte (bos) bir surumunu
# sys.modules'e yerlestiriyoruz -- Bio.PDB sadece bu isimlerin var olmasini
# kontrol ediyor, gercek islevini hic cagirmiyoruz. Bu satirlar Bio.PDB
# import edilmeden ONCE calismali.
fake_align = types.ModuleType("Bio.Align")
class _Dummy:
    pass
fake_align.Alignment = _Dummy
fake_align.MultipleSeqAlignment = _Dummy
fake_align.PairwiseAligner = _Dummy
sys.modules["Bio.Align"] = fake_align
 
from Bio.PDB import PDBParser, NeighborSearch
from Bio.SeqUtils import seq1
import warnings
warnings.filterwarnings("ignore")
 
# ---------------------------------------------------------------------------
# 1. AMINO ASIT OZELLIK TABLOLARI
# ---------------------------------------------------------------------------
 
# Kyte-Doolittle hidrofobiklik skalasi (1982) -- pozitif deger = daha hidrofobik
# Kaynak: Kyte J, Doolittle RF. J Mol Biol. 1982.
KYTE_DOOLITTLE = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}
 
# Amino asit hacimleri (A^3) -- Zamyatnin, 1972
RESIDUE_VOLUME = {
    "A":  88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G":  60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S":  89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}
 
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
 
 
# ---------------------------------------------------------------------------
# 2. GOMULULUK (BURIAL) HESABI
# ---------------------------------------------------------------------------
 
def compute_burial(chain, position, radius=10.0, normalization=30):
    """
    Bir pozisyonun ne kadar 'gomulu' oldugunu tahmin eder.
 
    YONTEM: "Contact Number" (temas sayisi) -- literaturde protein yapisi
    calismalarinda kullanilan, basit ama etkili bir yontem (orn. Echave ve
    calisma arkadaslarinin evrimsel kisitlamalar uzerine calismalari).
    SADECE CA (alfa-karbon) atomlarinin yogunluguna bakar -- yan zincirin
    yonune/uzunluguna bagli degildir, bu yuzden daha guvenilir bir olcumdur.
 
    Mantik: bir CA atominin belirli bir yaricap icinde ne kadar cok baska CA
    atomu varsa, o bolge o kadar sikica paketlenmis (gomulu) demektir.
    Yuzeydeki amino asitlerin etrafinda daha az CA atomu bulunur cunku bir
    tarafları cozucuye (suya) acik, "bos" bir alana bakar.
 
    normalization: kac komsu CA'nin 'tamamen gomulu' sayilacagini belirleyen
    kaba bir esik degeri (bu yapida deneyerek kalibre edildi, evrensel bir
    sabit degil -- farkli proteinlerde yeniden ayarlanmasi gerekebilir).
    """
    ca_atoms = [res["CA"] for res in chain if res.id[0] == " " and "CA" in res]
    ns = NeighborSearch(ca_atoms)
 
    target_ca = chain[position]["CA"]
    close = ns.search(target_ca.coord, radius)
    raw_count = len(close) - 1  # kendisini cikar (0 mesafede kendisiyle eslesir)
 
    burial = min(1.0, raw_count / normalization)  # 0 (yuzey) - 1 (gomulu) arasina sikistir
    return burial, raw_count
 
 
# ---------------------------------------------------------------------------
# 3. DDG_PROXY HESABI
# ---------------------------------------------------------------------------
 
def ddg_proxy(chain, position, wt_aa, mut_aa, w_hydro=0.05, w_vol=0.01):
    """
    Basit, seffaf DeltaDeltaG proxy skoru hesaplar.
    Pozitif = kararsizlastirici, Negatif = kararli artirici.
    """
    # Yapidaki gercek amino asidi dogrula
    actual_resname = chain[position].resname
    actual_aa = THREE_TO_ONE.get(actual_resname, "X")
    if actual_aa != wt_aa:
        print(f"UYARI: Pozisyon {position} icin beklenen '{wt_aa}', "
              f"yapida bulunan '{actual_aa}' ({actual_resname})")
 
    burial, raw_count = compute_burial(chain, position)
 
    delta_hydro = KYTE_DOOLITTLE[wt_aa] - KYTE_DOOLITTLE[mut_aa]  # kayip pozitifse hidrofobiklik azaliyor
    delta_vol = abs(RESIDUE_VOLUME[wt_aa] - RESIDUE_VOLUME[mut_aa])
 
    score = burial * (delta_hydro * w_hydro + delta_vol * w_vol)
 
    return {
        "position": position,
        "wt_aa": wt_aa,
        "mut_aa": mut_aa,
        "burial": round(burial, 3),
        "raw_neighbor_count": raw_count,
        "delta_hydrophobicity": round(delta_hydro, 2),
        "delta_volume": round(delta_vol, 2),
        "ddg_proxy": round(score, 3),
    }
 
 
def interpret(score):
    if score > 0.3:
        return "Kararsizlastirici (destabilizing) -- pozitif ddG, yapi zayifliyor beklentisi"
    elif score < -0.3:
        return "Kararli artirici (stabilizing) -- negatif ddG"
    else:
        return "Onemli bir etki beklenmiyor (notr)"
 
 
# ---------------------------------------------------------------------------
# 4. CALISTIR: I56T mutasyonu icin
# ---------------------------------------------------------------------------
 
if __name__ == "__main__":
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1REX", "data/1REX.pdb")
    chain_a = structure[0]["A"]
 
    print("=" * 70)
    print("DDG_PROXY HESABI: I56T mutasyonu (Izolosin -> Treonin)")
    print("=" * 70)
 
    result = ddg_proxy(chain_a, position=56, wt_aa="I", mut_aa="T")
 
    print(f"\nPozisyon           : {result['position']}")
    print(f"Wild-type amino asit: {result['wt_aa']} (Izolosin)")
    print(f"Mutant amino asit    : {result['mut_aa']} (Treonin)")
    print(f"Gomululuk skoru      : {result['burial']}  (0=yuzey, 1=tam gomulu)")
    print(f"  (ham komsu atom sayisi: {result['raw_neighbor_count']})")
    print(f"Hidrofobiklik kaybi  : {result['delta_hydrophobicity']}  (Kyte-Doolittle birimi)")
    print(f"Hacim farki          : {result['delta_volume']} A^3")
    print(f"\n>>> ddG_proxy skoru: {result['ddg_proxy']}")
    print(f">>> Yorum: {interpret(result['ddg_proxy'])}")
 
    # --- Karsilastirma: ayni analizi D67H icin de yapalim (referans mutasyon) ---
    print("\n" + "-" * 70)
    print("KARSILASTIRMA: D67H mutasyonu (alternatif, literaturde bilinen)")
    print("-" * 70)
    result2 = ddg_proxy(chain_a, position=67, wt_aa="D", mut_aa="H")
    print(f"\nPozisyon           : {result2['position']}")
    print(f"Gomululuk skoru      : {result2['burial']}")
    print(f"Hidrofobiklik kaybi  : {result2['delta_hydrophobicity']}")
    print(f"Hacim farki          : {result2['delta_volume']}")
    print(f"\n>>> ddG_proxy skoru: {result2['ddg_proxy']}")
    print(f">>> Yorum: {interpret(result2['ddg_proxy'])}")
 
    # --- Negatif kontrol: proteinin yuzeyinde, onemsiz beklenen bir mutasyon ---
    # Pozisyon 1 (Lys, dizinin en basi) genelde cok esnek/yuzeydedir.
    print("\n" + "-" * 70)
    print("NEGATIF KONTROL: K1A mutasyonu (yuzeyde, onemsiz beklenir)")
    print("-" * 70)
    result3 = ddg_proxy(chain_a, position=1, wt_aa="K", mut_aa="A")
    print(f"\nGomululuk skoru: {result3['burial']}")
    print(f">>> ddG_proxy skoru: {result3['ddg_proxy']}")
    print(f">>> Yorum: {interpret(result3['ddg_proxy'])}")
 