from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial import cKDTree

# Dictionary of van der Waals radii (in angstroms)
VDW_RADII = {
    'H': 1.10, 'B':1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'Si':2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.83, 'I': 1.98
}

def get_vdw_radius(symbol):
    return VDW_RADII.get(symbol, 1.5)  # Use default if not found

def generate_conformer(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def get_atom_coords_and_radii(mol):
    conf = mol.GetConformer()
    coords = []
    radii = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])
        radii.append(get_vdw_radius(atom.GetSymbol()))
    return np.array(coords), np.array(radii)

def compute_vdw_volume(coords, radii, grid_spacing=0.5):
    min_corner = coords.min(axis=0) - radii.max()
    max_corner = coords.max(axis=0) + radii.max()
    
    x, y, z = [np.arange(min_corner[i], max_corner[i], grid_spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(x, y, z, indexing='ij'), -1).reshape(-1, 3)
    
    tree = cKDTree(coords)
    nearby_atoms = tree.query_ball_point(grid, r=radii.max())

    occupied = 0
    for i, atom_indices in enumerate(nearby_atoms):
        point = grid[i]
        for j in atom_indices:
            if np.linalg.norm(point - coords[j]) <= radii[j]:
                occupied += 1
                break

    voxel_volume = grid_spacing ** 3
    return occupied * voxel_volume

def get_vdw_volume_from_smiles(smiles, grid_spacing=0.1):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = generate_conformer(mol)
    coords, radii = get_atom_coords_and_radii(mol)
    vol = compute_vdw_volume(coords, radii, grid_spacing)
    radius = ((3*vol)/(4*np.pi))**(1/3)
    return radius

if __name__ == "__main__":
    #molecule must be in SMILES Format
    moleculelist = [
    'F[B-]([N+]1=CC=CC1=C2)(F)N3C2=CC=C3',
    ]
    for name in moleculelist:
        print(get_vdw_volume_from_smiles(name))