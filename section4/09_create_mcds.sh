# Define your variables as before
ALLC_TABLE="/ceph/MethDev/JW240627--at-snmCT_with_TE/mCT_with_TE/stats/AllcPaths.tsv"
CHROM_SIZES="/ceph/MethDev/pbio/andy/JW/section_4/data/TAIR10_numeric.genome"
DMW_BED="/ceph/MethDev/pbio/andy/JW/section_4/data/chunks/chunks_CG_minfilt.fixed.bed"
GBM_BED="/ceph/MethDev/pbio/andy/JW/section_4/data/chunks/chunks_CG_regular_minfilt.fixed.bed"

# Note: generate-dataset prefers an output directory path, typically ending in .mcds
OUT_PATH="/ceph/MethDev/pbio/kay/section4/CG_chunks_mcds.mcds"

# Run the modern command
allcools generate-dataset \
  --allc_table "$ALLC_TABLE" \
  --output_path "$OUT_PATH" \
  --chrom_size_path "$CHROM_SIZES" \
  --obs_dim cell \
  --cpu 30 \
  --regions DMW "$DMW_BED" \
  --regions GbM "$GBM_BED" \
  --quantifiers DMW count CGN,CHN,CHG,CHH \
  --quantifiers GbM count CGN,CHN,CHG,CHH
