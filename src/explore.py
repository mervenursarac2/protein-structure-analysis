from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

parser = PDBParser(QUIET=True)
structure = parser.get_structure("1REX", "../data/1REX.pdb")

model = structure[0]
print("Zincir sayısı:", len(model))

for chain in model:
    print("Zincir ID:", chain.id) 


chain_a = model["A"]
sequence = ""

for residue in chain_a:
    if residue.id[0] == " ":  # su moleküllerini (HETATM) atla, sadece amino asitleri al
        sequence += seq1(residue.resname)

print("Toplam amino asit sayısı:", len(sequence))
print("Dizi:", sequence)

if sequence[55] == "I":
    print("Doğrulandı: pozisyon 56 İzolösin (I)")
else:
    print("Beklenmedik durum, bulunan:", sequence[55])
