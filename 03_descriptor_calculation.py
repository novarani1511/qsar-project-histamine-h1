import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np

def calculate_descriptors():
    df = pd.read_csv('h1_antagonists_clean.csv')
    print(f"Total senyawa yang akan diproses: {len(df)}")
    
    # Menyiapkan list semua deskriptor 2D di RDKit (sekitar 208 deskriptor)
    desc_names = [d[0] for d in Descriptors._descList]
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
    
    descriptors_list = []
    valid_smiles = []
    valid_pIC50 = []
    
    for i, row in df.iterrows():
        smiles = row['canonical_smiles']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is not None:
            # --- 1. Pemodelan Molekuler & Optimasi Geometri 3D ---
            # Menambahkan Hidrogen untuk optimasi 3D yang akurat
            mol_3d = Chem.AddHs(mol)
            # Generate konformer 3D awal
            embed_status = AllChem.EmbedMolecule(mol_3d, randomSeed=42, useRandomCoords=True)
            if embed_status == 0:
                # Minimisasi energi (Optimasi Geometri) dengan MMFF94 force field
                try:
                    AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
                except Exception as e:
                    pass # Abaikan jika gagal konvergensi dalam 200 iterasi
            
            # --- 2. Perhitungan Deskriptor ---
            # RDKit descriptors umumnya lebih konsisten jika dihitung tanpa explicit Hs (tergantung spesifikasi deskriptor)
            mol_for_desc = Chem.RemoveHs(mol_3d)
            
            try:
                desc_values = calculator.CalcDescriptors(mol_for_desc)
                descriptors_list.append(desc_values)
                valid_smiles.append(smiles)
                valid_pIC50.append(row['pIC50'])
            except Exception as e:
                print(f"Gagal menghitung deskriptor untuk SMILES index {i}")
                
        if (i+1) % 50 == 0:
            print(f"{i+1} senyawa telah dioptimasi dan dihitung deskriptornya...")

    print("Selesai menghitung deskriptor.")
    
    # Gabungkan menjadi DataFrame
    df_desc = pd.DataFrame(descriptors_list, columns=desc_names)
    df_desc['canonical_smiles'] = valid_smiles
    df_desc['pIC50'] = valid_pIC50
    
    # Pindahkan SMILES dan pIC50 ke bagian depan DataFrame
    cols = ['canonical_smiles', 'pIC50'] + desc_names
    df_desc = df_desc[cols]
    
    # Simpan
    output_file = 'h1_antagonists_descriptors.csv'
    df_desc.to_csv(output_file, index=False)
    print(f"Data deskriptor berhasil disimpan ke {output_file} dengan dimensi {df_desc.shape}")

if __name__ == "__main__":
    calculate_descriptors()
