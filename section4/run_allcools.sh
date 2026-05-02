ALLC_TABLE="/ceph/MethDev/JW240627--at-snmCT_with_TE/mCT_with_TE/stats/AllcPaths.tsv"
GENOME="/gale/netapp/home/kqu/allcools_work/data/genome.sizes"
DMW="/gale/netapp/home/kqu/allcools_work/data/DMW_clean_3col.bed"
GBM="/gale/netapp/home/kqu/allcools_work/data/GbM_clean_3col.bed"
OUT_PATH="/gale/netapp/home/kqu/allcools_work/CG_chunks.mcds"
rm -rf "$OUT_PATH"
allcools generate-dataset \
  --allc_table "$ALLC_TABLE" \
  --output_path "$OUT_PATH" \
  --chrom_size_path "$GENOME" \
  --obs_dim cell \
  --cpu 10 \
  --regions DMW "$DMW" \
  --regions GbM "$GBM" \
  --quantifiers DMW count CGN,CHN,CHG,CHH \
  --quantifiers GbM count CGN,CHN,CHG,CHH
