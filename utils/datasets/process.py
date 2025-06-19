from Bio import SeqIO
from Bio.PDB import PDBParser, Polypeptide
from Bio.PDB.Polypeptide import PPBuilder
from Bio import PDB
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import warnings
warnings.filterwarnings("ignore")
# 提取蛋白、多肽序列
# 定义三字母到一字母的映射字典
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequences(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    sequences = []
    for model in structure:
        for chain in model:
            seq = ''
            for residue in chain:
                # 获取三字母代码并转换为一字母代码
                resname = residue.resname
                if resname in three_to_one:
                    seq += three_to_one[resname]
                else:
                    seq += 'X'  # 未知氨基酸标记为 'X'
            sequences.append(seq)
    return sequences
# def extract_sequences(pdb_file):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure('protein', pdb_file)
#     sequences = []

#     for model in structure:
#         for chain in model:
#             # 提取链的氨基酸序列
#             polypeptides = Polypeptide.PPBuilder().build_peptides(chain)
#             for poly_index, poly_index in enumerate(polypeptides):
#                 # 转换为一字母表示法
#                 seq = poly_index.get_sequence()
#                 sequences.append(str(seq))

#     return sequences
# 提取蛋白的基本3D结构信息
def extract_pdb_features(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    features = []

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                resseq = residue.get_id()[1]  # 残基序列号
                resname = residue.get_resname()  # 氨基酸类型
                atom_coords = []
                for atom in residue:
                    atom_coords.append(atom.get_coord())  # 原子坐标
                features.append({
                    'chain_id': chain_id,
                    'resseq': resseq,
                    'res_nb': residue.get_id(),
                    'aa': resname,
                    'pos_heavyatom': atom_coords,
                })
    return features
# 提取B因子、icode
def extract_bfactors(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    features = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    features.append({
                        'icode': residue.get_id()[2],  # 插入代码
                        'bfactor_heavyatom': atom.get_bfactor()  # B因子
                    })
    return features
# 计算蛋白质的主链二面角(phi, psi, omega)
def extract_dihedral_angles(pdb_file):
    u = mda.Universe(pdb_file)
    rama = Ramachandran(u.select_atoms('backbone'))
    rama.run()

    dihedral_angles = {'phi': [], 'psi': [], 'omega': []}

    # 提取 phi 和 psi 角度
    for ts in rama.angles:
        phi_psi = ts
        dihedral_angles['phi'].append(phi_psi[0])
        dihedral_angles['psi'].append(phi_psi[1])

    # Omega 角度需要额外计算
    backbone = u.select_atoms('backbone')
    for i in range(1, len(backbone.residues) - 1):
        prev_residue = backbone.residues[i - 1]
        current_residue = backbone.residues[i]
        next_residue = backbone.residues[i + 1]
        omega = calculate_omega(prev_residue, current_residue, next_residue)
        dihedral_angles['omega'].append(omega)

    return dihedral_angles

def get_coords(residue, atom_name):
    """获取特定原子的坐标。"""
    atom = residue.atoms.select_atoms(f'name {atom_name}')
    if len(atom) > 0:
        return atom.positions[0]
    else:
        return None

def calculate_omega(residue1, residue2, residue3):
    """计算 omega 角。"""
    atoms = ['N', 'CA', 'C']
    coords = []

    for atom_name in atoms:
        for residue in [residue1, residue2, residue3]:
            coord = get_coords(residue, atom_name)
            if coord is not None:
                coords.append(coord)
    
    if len(coords) < 4:
        raise ValueError("没有足够的坐标来计算 omega 角")

    coords = np.array(coords)

    def dihedral(p1, p2, p3, p4):
        """计算二面角。"""
        b0 = p2 - p1
        b1 = p3 - p2
        b2 = p4 - p3
        v = np.cross(b0, b1)
        w = np.cross(b1, b2)
        x = np.dot(v, w)
        y = np.dot(np.cross(b0, b1), b2) / np.linalg.norm(b1)
        return np.arctan2(x, y)
    
    omega = dihedral(coords[0], coords[1], coords[2], coords[3])
    return np.degrees(omega)

# 计算蛋白的侧链的扭转角度
def extract_chi_angles(pdb_file):
    u = mda.Universe(pdb_file)
    chi_angles = []
    
    # Ramachandran分析器用于提取二面角，包括Chi角度
    rama = Ramachandran(u.select_atoms('protein'))
    rama.run()
    
    for ts in rama.angles:
        chi_angles.append({
            'chi1': ts[0],
            'chi2': ts[1] if len(ts) > 1 else None,
            'chi_complete': all([v is not None for v in ts])  # 检查Chi角度是否完整
        })
    
    return chi_angles
# 3D结构特征封装
def integrate_features(pdb_file):
    sequence_features = extract_pdb_features(pdb_file)
    bfactor_features = extract_bfactors(pdb_file)
    dihedral_angles = extract_dihedral_angles(pdb_file)
    chi_angles = extract_chi_angles(pdb_file)
    
    # 检查各个特征的长度
    print(f"Number of sequence_features: {len(sequence_features)}")
    print(f"Number of bfactor_features: {len(bfactor_features)}")
    print(f"Number of dihedral_angles (phi): {len(dihedral_angles['phi'])}")
    print(f"Number of chi_angles: {len(chi_angles)}")
    
    # 保证所有列表长度一致
    min_length = min(len(sequence_features), len(bfactor_features), len(dihedral_angles['phi']), len(chi_angles))
    
    integrated_features = []
    for i in range(min_length):
        feature_dict = {
            **sequence_features[i],
            **(bfactor_features[i] if i < len(bfactor_features) else {}),
            'phi': dihedral_angles['phi'][i] if i < len(dihedral_angles['phi']) else None,
            'psi': dihedral_angles['psi'][i] if i < len(dihedral_angles['psi']) else None,
            'omega': dihedral_angles['omega'][i] if i < len(dihedral_angles['omega']) else None,
            **(chi_angles[i] if i < len(chi_angles) else {}),
        }
        integrated_features.append(feature_dict)
    
    return integrated_features
# 提取多肽和蛋白的序列
peptide_seq = extract_sequences('/home/zht/programs/PP_generate/data/ppbench2024/1a0m_A/peptide.pdb')
receptor_seq = extract_sequences('/home/zht/programs/PP_generate/data/ppbench2024/1a0m_A/receptor.pdb')
print(f"Peptide Sequence: {peptide_seq}")
print(f"Receptor Sequence: {receptor_seq}")
# 3D结构特征封装
# peptide_integrated_features = integrate_features('/home/zht/programs/PP_generate/data/ppbench2024/1a0m_A/peptide.pdb')
# receptor_integrated_features = integrate_features('/home/zht/programs/PP_generate/data/ppbench2024/1a0m_A/receptor.pdb')
# 打印整合的特征
# for feature in peptide_integrated_features[:3]:
#     print(feature)






