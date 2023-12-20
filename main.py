from Bio.PDB import *
import sys

promienie = {"H": 1.2, "O": 1.52, "N": 1.55, "C": 1.7} #kąty sil van der waalsa

file = sys.argv[1]

pdb_parser = PDBParser()
structure = pdb_parser.get_structure("model", file)
bad_overlaps = 0
atoms = list(structure.get_atoms())
no_atoms = len(atoms)
odwiedzone = set()
neighbor_search = NeighborSearch(atoms)

for atom in atoms:
    res = atom.get_parent()
    for atom1 in neighbor_search.search(atom.coord, 4, level='A'):
        res1 = atom1.get_parent()

        atom_pair = frozenset([atom, atom1])

        #warunek sprawdzajacy czy para atomow zostala juz wzieta pod uwage
        if atom_pair in odwiedzone:
            continue
        odwiedzone.add(atom_pair)

        if abs(res.id[1] - res1.id[1]) < 2 or res == res1:
            continue

        distance = atom - atom1

        promien1 = promienie.get(atom.element, 1.5)
        promien2 = promienie.get(atom1.element, 1.5)


        #sprawdzanie czy promienie dzialania sil na siebie nachodza
        if (promien1 + promien2 - distance) >= 0.4:
            bad_overlaps += 1


#wzor do liczeniu clash score
clash_score = (1000 * bad_overlaps) / no_atoms

print("Clash Score: ", clash_score)