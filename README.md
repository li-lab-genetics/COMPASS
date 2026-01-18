# COMPASS: COding variant Mapping onto Protein-structure-Anchored hotSpotS
A computational pipeline for structural and functional interpretation of disease-associated coding variants.
## Description
COMPASS is a computational framework designed to integrate whole-genome sequencing (WGS) data with experimentally determined and AI-predicted protein structures to systematically link disease-associated coding variants with functional structural regions and therapeutic targets. By combining sequence-based variant selection and structure-based patch scanning, COMPASS identifies disease-relevant substructures (functional hotspots) that accelerate genetics-driven drug discovery.
## Workflow Overview
The framework comprises three modules: (i) a sequence module (COMPASS-Seq) that maps variants to transcript-specific amino acid sequences and consolidates association signals at the residue level; (ii) a structure module (COMPASS-Struct) that localizes disease-relevant hotspots in the protein; and (iii) an annotation module (COMPASS-Anno) that integrates functional annotations and automated literature mining to support interpretation.


![COMPASS_workflow](docs/workflow.jpg)
## Docker Image
A Docker image for COMPASS, which includes both R and Python environments with all COMPASS-related dependencies pre-installed, is available on Docker Hub. The image can be pulled using:
```
docker pull fengyn923/compass:0.1.0
```
## Examples
### COMPASS-Seq
```
Rscript coding_variants2amino_acids.R \
  --category=missense \
  --obj_nullmodel.file=obj.STAAR.UKB.Brain_Spine_Secondary.20250422.Rdata \
  --chr=6 \
  --gds.path=ukb.500k.wgs.chr6.pass.annotated.extend.gds \
  --gene_name=CRIP3 \
  --transcript=ENST00000274990 \
  --protein_sequence_selection_strategy=one-step
```
| Parameter                               | Example                                              | Description                                                                                                                                       |
| --------------------------------------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--category`                            | `missense`                                           | Type of coding variants to analyze. Options: `missense`, `disruptive missense`, `ptv`, `plof`.                                                    |
| `--obj_nullmodel.file`                  | `obj.STAAR.UKB.Brain_Spine_Secondary.20250422.Rdata` | The pre-computed null model object was generated for the target phenotype.                                                                        |
| `--chr`                                 | `6`                                                  | Chromosome number containing the target gene.                                                                                                     |
| `--gds.path`                            | `ukb.500k.wgs.chr6.pass.annotated.extend.gds`        | Path to the annotated GDS file containing WGS variants for the specified chromosome.                                                              |
| `--gene_name`                           | `CRIP3`                                              | Name of the target gene for analysis.                                                                                                             |
| `--transcript`     | `ENST00000333681`                                        | Ensembl transcript ID corresponding to the analyzed gene. The default setting is `standard`, which automatically selects the canonical transcript for analysis.|
| `--protein_sequence_selection_strategy` | `one-step`                                           | Strategy used to select representative amino acid sequences when multiple variants occur at the same position. Options: `one-step`, `stepwise`, `best-subset-selection`. |

### COMPASS-Struct
```
python patch_scanning.py \
  --i_pdb_file bcl2_chronic_lymphoid_leukemia.pdb \
  --i_variant_file BCL2_ENST00000333681_onestep_variant.txt \
  --radius 20 \
  --transcript ENST00000333681 \
  --gene_name BCL2 \
  --obj_nullmodel obj.STAAR.UKB.Chronic_Lymphoid_Leukemia.20250422.Rdata \
  --chr 18 \
  --category missense \
  --gds_path ukb.500k.wgs.chr18.pass.annotated.extend.gds \
  --num_patch 10
```

| Parameter          | Example                                                  | Description                                                                                                     |
| ------------------ | -------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `--i_pdb_file`     | `bcl2_chronic_lymphoid_leukemia.pdb`                     | Input PDB file containing the 3D structure of the target protein.                                               |
| `--i_variant_file` | `BCL2_ENST00000333681_onestep_variant.txt`               | Variant file listing disease-associated variants (output from the sequence-based selection step).               |
| `--radius`         | `20`                                                     | Radius (in Å) used to define a 3D patch centered on each residue’s Cα atom.                                     |
| `--transcript`     | `ENST00000333681`                                        | Ensembl transcript ID corresponding to the analyzed gene. The default setting is `standard`, which automatically selects the canonical transcript for analysis.|
| `--gene_name`      | `BCL2`                                                   | Name of the target gene for analysis.                                                                           |
| `--obj_nullmodel`  | `obj.STAAR.UKB.Chronic_Lymphoid_Leukemia.20250422.Rdata` | The pre-computed null model object was generated for the target phenotype.                                      |
| `--chr`            | `18`                                                     | Chromosome number containing the target gene.                                                                   |
| `--category`       | `missense`                                               | Variant category to analysis in the patch scanning.                                                             |
| `--gds_path`       | `ukb.500k.wgs.chr18.pass.annotated.extend.gds`           | Path to the annotated GDS file containing WGS variants for the specified chromosome.                            |
| `--num_patch`      | `10`                                                     | Specifies the number of top-ranked structural hotspots (patches) with the lowest P-values to be reported. A value of 999 outputs all identified hotspots. |

### COMPASS-Anno
COMPASS-Anno is used to provide functional interpretation for disease-associated structural hotspots identified by COMPASS. The interpretation is achieved by integrating functional context and literature-based evidence, with assistance from large language models, including ChatGPT and DeepSeek, and is accessible through our website (http://compassbio.net).



## Authors
Yannuo Feng<sup>&ast;</sup>, Yihao Peng<sup>&ast;</sup>, Shijie Fan<sup>&ast;</sup>, Shijia Bian, Jingjing Gong, Qikai Jiang, Chang Lu<sup>#</sup>, Xihao Li<sup>#</sup>, Zilin Li<sup>#</sup>. **Prioritizing druggable targets by mapping human disease-associated coding variants onto protein structures with COMPASS**.



Contact: luc816@nenu.edu.cn, xihaoli@unc.edu, lizl@nenu.edu.cn

