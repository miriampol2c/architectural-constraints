import os
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, PDBIO

def cif_to_pdb(cif_path, pdb_path):
    """Convert an mmCIF file to PDB using Biopython."""
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    try:
        structure = parser.get_structure(os.path.basename(cif_path), cif_path)
        io.set_structure(structure)
        io.save(pdb_path)
        print(f"Converted: {os.path.basename(cif_path)} → {os.path.basename(pdb_path)}")
    except Exception as e:
        print(f"Failed: {os.path.basename(cif_path)} → {e}")

def convert_all_cifs(input_dir, output_dir, num_workers=4):
    """Convert all CIF files in input_dir to PDB format in output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    cif_files = [f for f in os.listdir(input_dir) if f.endswith(".cif")]
    if not cif_files:
        print(f"No .cif files found in {input_dir}")
        return

    tasks = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for cif_file in cif_files:
            cif_path = os.path.join(input_dir, cif_file)
            pdb_path = os.path.join(output_dir, os.path.splitext(cif_file)[0] + ".pdb")
            tasks.append(executor.submit(cif_to_pdb, cif_path, pdb_path))

        for task in tasks:
            task.result()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert all CIF files in a directory to PDB format.")
    parser.add_argument("input_dir", help="Directory containing CIF files")
    parser.add_argument("output_dir", help="Directory to save converted PDB files")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (default: 4)")
    args = parser.parse_args()

    convert_all_cifs(args.input_dir, args.output_dir, args.workers)

