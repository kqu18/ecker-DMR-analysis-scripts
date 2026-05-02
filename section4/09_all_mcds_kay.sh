#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=24:00:00
#$ -l mem_free=32G
#$ -pe smp 30
#$ -N mcds_chunks
#$ -o /ceph/MethDev/pbio/andy/JW/section_4/logs/mcds_chunks.out
#$ -e /ceph/MethDev/pbio/andy/JW/section_4/logs/mcds_chunks.err

set -euo pipefail

SCRIPT_DIR="/ceph/MethDev/pbio/andy/JW/section_4"
KAY_DIR="${SCRIPT_DIR}/kay"

mkdir -p "${SCRIPT_DIR}/logs"
cd "${KAY_DIR}"

echo "Running ${KAY_DIR}/09_all_mcds_kay.sh"
echo "Job started at: $(date)"

ALLC_TABLE="/ceph/MethDev/JW240627--at-snmCT_with_TE/mCT_with_TE/stats/AllcPaths.tsv"
CHROM_SIZES="/gale/raidix/rdx-7/tnobori/tools/YAP/reference_files/Arabidopsis_thaliana.TAIR10.dna.toplevel_chrL_appended_sizes.genome"

DMW_BED="${SCRIPT_DIR}/data/chunks/chunks_CG_minfilt.fixed.bed.gz"
GBM_BED="${SCRIPT_DIR}/data/chunks/chunks_CG_regular_minfilt.fixed.bed.gz"

OUT_PREFIX="${SCRIPT_DIR}/CG_chunks_mcds_v3"

rm -rf "${OUT_PREFIX}.tmp_dir" "${OUT_PREFIX}_0.tmp_dir" "${OUT_PREFIX}_1.tmp_dir"


source /gale/netapp/home/kqu/miniconda3/etc/profile.d/conda.sh


echo "CONDA_PREFIX=$CONDA_PREFIX"
which allcools
python -c "import pandas; print(pandas.__version__)"

conda activate allcools_kay
which allcools
python -c "import pandas; print(pandas.__version__)"


allcools generate-mcds \
  --allc_table "${ALLC_TABLE}" \
  --output_prefix "${OUT_PREFIX}" \
  --chrom_size_path "${CHROM_SIZES}" \
  --mc_contexts CGN CHN CHG CHH \
  --region_bed_paths "${DMW_BED}" "${GBM_BED}" \
  --region_bed_names DMW GbM \
  --cpu 30 \
  --max_per_mcds 100000

echo "Job finished at: $(date)"
