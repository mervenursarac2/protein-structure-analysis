from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
structure = parser.get_structure("1REX", "../data/1REX.pdb")
model = structure[0]
chain_a = model["A"]

# I56'nın kendisini bul
target_residue = chain_a[56]
print("Hedef amino asit:", target_residue.resname, target_residue.id[1])

# Yapıdaki TÜM atomları topla (komşu aramasi icin gerekli)
all_atoms = [atom for atom in chain_a.get_atoms()]
ns = NeighborSearch(all_atoms)

# I56'nın yan zincirindeki atomlara (CB, CG1, CG2, CD1 -- Izolösin'in "kolları") bak
# ve her birinin 5 Angstrom yakınındaki atomları bul
nearby_residues = set()

for atom in target_residue:
    if atom.get_name() not in ("N", "CA", "C", "O"):  # sadece yan zincir atomları
        close_atoms = ns.search(atom.coord, 5.0)  # 5 Angstrom yarıçapında ara
        for close_atom in close_atoms:
            parent_residue = close_atom.get_parent()
            if parent_residue.id[1] != 56:  # kendisini sayma
                nearby_residues.add((parent_residue.id[1], parent_residue.resname))

# Sonuçları pozisyona göre sırala ve yazdır
print("\nI56'nın 5 Angstrom yakınındaki amino asitler:")
for pos, resname in sorted(nearby_residues):
    print(f"  {pos}: {resname}")