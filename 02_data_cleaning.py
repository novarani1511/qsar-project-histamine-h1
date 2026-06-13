import pandas as pd
import numpy as np

def clean_data():
    input_file = 'h1_antagonists_raw.csv'
    df = pd.read_csv(input_file)
    print(f"Data awal: {len(df)} baris")
    
    # Filter data dengan nilai IC50 dan SMILES yang valid
    df = df[df['standard_value'].notna()]
    df = df[df['canonical_smiles'].notna()]
    print(f"Setelah filter missing value: {len(df)} baris")
    
    # Pastikan unit dalam nM
    df = df[df['standard_units'] == 'nM']
    print(f"Setelah filter unit nM: {len(df)} baris")
    
    # Ubah ke numerik
    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df.dropna(subset=['standard_value'])
    
    # Filter nilai 0 atau negatif (karena akan di-log)
    df = df[df['standard_value'] > 0]
    
    # Hitung pIC50
    # pIC50 = -log10(IC50 in M) = -log10(IC50 in nM * 1e-9) = 9 - log10(IC50 in nM)
    df['pIC50'] = 9 - np.log10(df['standard_value'])
    
    # Bersihkan duplikat berdasarkan SMILES, ambil rata-rata pIC50
    df_clean = df.groupby('canonical_smiles').agg({'pIC50': 'mean'}).reset_index()
    print(f"Setelah menghapus duplikat: {len(df_clean)} senyawa unik")
    
    output_file = 'h1_antagonists_clean.csv'
    df_clean.to_csv(output_file, index=False)
    print(f"Data bersih berhasil disimpan ke {output_file}")

if __name__ == "__main__":
    clean_data()
