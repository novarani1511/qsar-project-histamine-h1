import subprocess
import os
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def prepare_and_run_vina():
    print("--- 1. Persiapan Ligan ---")
    # Ambil senyawa paling aktif dari data yang kita bersihkan
    df = pd.read_csv('h1_antagonists_clean.csv')
    df = df.sort_values(by='pIC50', ascending=False)
    
    # Ambil senyawa dengan pIC50 tertinggi sebagai contoh ligan yang akan di-docking
    best_smiles = df.iloc[0]['canonical_smiles']
    best_pic50 = df.iloc[0]['pIC50']
    print(f"Senyawa aktif terbaik (pIC50 = {best_pic50:.2f}): {best_smiles}")
    
    # Generate 3D dari SMILES
    mol = Chem.MolFromSmiles(best_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Konversi ke PDBQT menggunakan library python 'meeko'
    try:
        from meeko import MoleculePreparation
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        ligand_pdbqt_string = preparator.write_pdbqt_string()
        with open("ligand.pdbqt", "w") as f:
            f.write(ligand_pdbqt_string)
        print("Ligan berhasil disiapkan dan disimpan sebagai 'ligand.pdbqt'")
    except ImportError:
        print("Modul 'meeko' tidak ditemukan. Harap jalankan: pip install meeko")
        return

    print("\n--- 2. Persiapan Reseptor ---")
    pdb_file = "3rze.pdb"
    if not os.path.exists(pdb_file):
        print(f"Mengunduh struktur kristal H1 (PDB ID: 3RZE)...")
        urllib.request.urlretrieve("https://files.rcsb.org/download/3RZE.pdb", pdb_file)
        print(f"Struktur {pdb_file} berhasil diunduh.")

    receptor_pdbqt = "receptor.pdbqt"
    if not os.path.exists(receptor_pdbqt):
        print("\n[PERHATIAN]: File 'receptor.pdbqt' tidak ditemukan!")
        print("Anda harus melakukan preparasi protein secara manual menggunakan AutoDockTools:")
        print("  1. Buka 3rze.pdb di AutoDockTools.")
        print("  2. Hapus ligan bawaan, molekul air, dan rantai yang tidak diperlukan.")
        print("  3. Tambahkan atom hidrogen polar dan muatan Kollman.")
        print("  4. Simpan sebagai 'receptor.pdbqt' di folder ini.")
        print("Skrip akan berhenti di sini. Silakan jalankan kembali setelah 'receptor.pdbqt' tersedia.")
        return

    print("\n--- 3. Menjalankan AutoDock Vina ---")
    if not os.path.exists("vina.exe"):
        print("\n[PERHATIAN]: File 'vina.exe' tidak ditemukan!")
        print("Silakan unduh AutoDock Vina untuk Windows (vina.exe) dari https://github.com/ccsb-scripps/AutoDock-Vina/releases")
        print("dan letakkan file vina.exe di folder yang sama dengan skrip ini.")
        return

    # Menjalankan Vina melalui Python subprocess
    print("Mengeksekusi Vina...")
    # NOTE: Koordinat Grid Box di bawah ini hanyalah ilustrasi.
    # Anda harus menggunakan AutoDockTools untuk menemukan koordinat (center) dari ligan asli (Doxepin)
    cmd = [
        "vina.exe",
        "--receptor", receptor_pdbqt,
        "--ligand", "ligand.pdbqt",
        "--center_x", "15.0", 
        "--center_y", "20.0",
        "--center_z", "45.0",
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--out", "docking_result.pdbqt",
        "--log", "docking_log.txt"
    ]
    
    process = subprocess.run(cmd, capture_output=True, text=True)
    
    if process.returncode == 0:
        print("\nProses Docking Selesai!")
        print("Hasil tersimpan di 'docking_result.pdbqt'.")
        print("Log tersimpan di 'docking_log.txt'.")
    else:
        print(f"\nTerjadi kesalahan saat menjalankan Vina:\n{process.stderr}")

if __name__ == "__main__":
    prepare_and_run_vina()
