from Bio.PDB import PDBParser
from scipy.spatial import KDTree
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import argparse
import subprocess
import io

parser = argparse.ArgumentParser(description="Process gene variant CSV to transcript-based files.")
parser.add_argument('--i_pdb_file', required=True, help='Path to input PDB file.')
parser.add_argument('--i_variant_file', required=True, help='Path to input variant file.')
parser.add_argument('--radius', type=int, default=15, help='Radius for neighbor search.')
parser.add_argument('--gene_name', required=True, help='Gene name.')
parser.add_argument('--transcript', required=True, help='Transcript ID for filtering variants.')
parser.add_argument('--obj_nullmodel', required=True, help='Path to null_model file.')
parser.add_argument('--chr', type=int,  help='Chromosome number.')
parser.add_argument('--num_patch', default=10, help='Path to neighbor file.')
parser.add_argument('--category', required=True, help='Variant category.')
parser.add_argument('--gds_path', required=True, help='Path to the GDS file.')
parser.add_argument('--output_dir', help='Directory to save output files.')
args = parser.parse_args()

print(args)

# Read the CSV file
pdb_file = args.i_pdb_file
variant_file = args.i_variant_file
gene_name = args.gene_name
radius = args.radius
transcript = args.transcript
num_patch = args.num_patch

parser = PDBParser()
structure = parser.get_structure("protein", pdb_file)
ca_coords = []
residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            if "CA" in residue:
                ca_coords.append(residue["CA"].get_coord())
                residues.append((residue.get_resname(), residue.get_id()[1]))
ca_coords = np.array(ca_coords)

kdtree = KDTree(ca_coords)

results = []
neighbor_ids = defaultdict(list)

for i, coord in enumerate(ca_coords):
    neighbors_idx = kdtree.query_ball_point(coord, radius)
    for idx in neighbors_idx:
        if idx != i:
            dist = np.linalg.norm(ca_coords[i] - ca_coords[idx])
            results.append([
                residues[i][0],
                residues[i][1],
                residues[idx][0],
                residues[idx][1], 
                dist
            ])
            neighbor_ids[residues[i][1]].append(residues[idx][1])


df = pd.read_csv(variant_file, sep = "\t")
exonic_info = df['GENCODE.EXONIC.Info']

aa_change = {} 
for idx, info in enumerate(exonic_info):
    variants = info.split(',')
    for variant in variants:
        if transcript in variant:
            match = re.search(r'p\.([A-Z](\d+)[A-Z])', variant)
            if match:
                match = re.match(r'([A-Z])(\d+)([A-Z])', match.group(1))
                position = int(match.group(2))
                aa_change[position] = idx

output_lines = []
data = []
for central_id, neighbors in neighbor_ids.items():
    all_positions = []
    all_positions.append(central_id)
    all_positions.extend(neighbors)
    variant_indices = [str(aa_change[pos]) for pos in all_positions if pos in aa_change]
    variant_indices_plus_one = []
    for idx in variant_indices:
        variant_indices_plus_one.append(str(int(idx) + 1))
    variant_str = ', '.join(variant_indices_plus_one) if variant_indices_plus_one else ''
    neighbors_str = ', '.join(map(str, neighbors))
    data.append({
        'central_id': central_id,
        'neighbors': neighbors_str,
        'variant_indices': variant_str
    })

neighbor_with_variant = pd.DataFrame(data)
neighbor_tsv = neighbor_with_variant.to_csv(index=False, sep="\t")

rscript_path = "Rscript"
script_path = "patch_scanning_p.R"
r_args = [
    rscript_path,
    script_path,
    f"--obj_nullmodel.file={args.obj_nullmodel}",
    f"--variant_file={args.i_variant_file}",
    f"--chr={args.chr}",
    f"--gene_name={args.gene_name}",
    f"--category={args.category}",
    f"--gds.path={args.gds_path}" 
]

# run R script
try:
    result = subprocess.run(r_args, input=neighbor_tsv, capture_output=True, text=True, check=True)
    print("R script output:")
    print(result.stdout)
    if result.stderr:
        print("R script errors (STDERR):")
        print(result.stderr)
    lines = result.stdout.splitlines()
    data_lines = [line for line in lines if len(line.split("\t")) == 8]
    data_str = "\n".join(data_lines)
    patch_df = pd.read_csv(io.StringIO(data_str), sep="\t")

except subprocess.CalledProcessError as e:
    print("❌ R script execution failed!")

    print("\n---- STDOUT ----")
    print(e.stdout or "[No stdout]")

    print("\n---- STDERR ----")
    print(e.stderr or "[No stderr]")

    print("\n---- Exception Details ----")
    print(f"Return Code: {e.returncode}")
    print(f"Command: {e.cmd}")


last_col = patch_df.columns[-1]

patch_df_sorted = patch_df.sort_values(by=last_col, ascending=True)

num_patch = int(num_patch)
if num_patch == 999:
    top_patches = patch_df_sorted
else:
    top_patches = patch_df_sorted.head(num_patch)


third_col_name = top_patches.columns[2]

result_rows = []

for idx, row in top_patches.iterrows():

    values = [int(v.strip()) for v in str(row[third_col_name]).split(',')]
    values = [value - 1 for value in values]

    selected_rows = df.iloc[[index for index in values]]
    exonic_info = selected_rows['GENCODE.EXONIC.Info']

    aa_change = []
    for info in exonic_info:
        variants = info.split(',')
        for variant in variants:
            if transcript in variant:
                match = re.search(r'p\.([A-Z]\d+[A-Z])', variant)
                aa_change.append(match.group(1))

    new_row = row.to_dict()
    new_row['AA_Change'] = ','.join(aa_change)
    result_rows.append(new_row) 

result_df = pd.DataFrame(result_rows)
output_file_path = f'{gene_name}_{transcript}_patch.csv'
result_df.to_csv(output_file_path, index=False)



    

