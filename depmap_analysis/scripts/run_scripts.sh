# Set the following environment variables:
# OUTDIR : the output directory
# DATA_DIR : where the inpup depmap data is
# SIF_DIR : where the DB dumps are
# Remember to check the mu sigma values for any new depmap dataset

ll=1; ul=2;
echo "Running $ll-$ul SD"
echo "Output goes to ${OUTDIR}/${ll}_${ul}sd/"
python depmap_script.py -cf "${DATA_DIR}/19Q4/Achilles_gene_effect.csv" -cc "${OUTDIR}/_crispr_all_correlations.h5"\
  -rf "${DATA_DIR}/demeter/D2_combined_gene_dep_scores.csv" -rc "${OUTDIR}/_rnai_all_correlations.h5"\
  -b "${SIF_DIR}/belief_dict.pkl" --filter-type None --no-web-files\
  --explained-set /home/klas/repos/depmap_analysis/input_data/depmap/gene_sets/mito_function_depmap.csv gene\
  -lw "${SIF_DIR}/db_dump_df_lite.csv" -cstats 0.003076 0.056813 -rstats 0.006854 0.077614 -crange $ll $ul -rrange $ll $ul\
  -o "${OUTDIR}/${ll}_${ul}sd/" > "${OUTDIR}/${ll}_${ul}sd/stout.log"
echo "Finished running $ll-$ul SD"

ll=2; ul=3;
echo "Running $ll-$ul SD"
echo "Output goes to ${OUTDIR}/${ll}_${ul}sd/"
python depmap_script.py -cf "${DATA_DIR}/19Q4/Achilles_gene_effect.csv" -cc "${OUTDIR}/_crispr_all_correlations.h5"\
  -rf "${DATA_DIR}/demeter/D2_combined_gene_dep_scores.csv" -rc "${OUTDIR}/_rnai_all_correlations.h5"\
  -b "${SIF_DIR}/belief_dict.pkl" --filter-type None --no-web-files\
  --explained-set /home/klas/repos/depmap_analysis/input_data/depmap/gene_sets/mito_function_depmap.csv gene\
  -lw "${SIF_DIR}/db_dump_df_lite.csv" -cstats 0.003076 0.056813 -rstats 0.006854 0.077614 -crange $ll $ul -rrange $ll $ul\
  -o "${OUTDIR}/${ll}_${ul}sd/" > "${OUTDIR}/${ll}_${ul}sd/stout.log"
echo "Finished running $ll-$ul SD"

ll=3; ul=4;
echo "Running $ll-$ul SD"
echo "Output goes to ${OUTDIR}/${ll}_${ul}sd/"
python depmap_script.py -cf "${DATA_DIR}/19Q4/Achilles_gene_effect.csv" -cc "${OUTDIR}/_crispr_all_correlations.h5"\
  -rf "${DATA_DIR}/demeter/D2_combined_gene_dep_scores.csv" -rc "${OUTDIR}/_rnai_all_correlations.h5"\
  -b "${SIF_DIR}/belief_dict.pkl" --filter-type None --no-web-files\
  --explained-set /home/klas/repos/depmap_analysis/input_data/depmap/gene_sets/mito_function_depmap.csv gene\
  -lw "${SIF_DIR}/db_dump_df_lite.csv" -cstats 0.003076 0.056813 -rstats 0.006854 0.077614 -crange $ll $ul -rrange $ll $ul\
  -o "${OUTDIR}/${ll}_${ul}sd/" > "${OUTDIR}/${ll}_${ul}sd/stout.log"
echo "Finished running $ll-$ul SD"

ll=4; ul=5;
echo "Running $ll-$ul SD"
echo "Output goes to ${OUTDIR}/${ll}_${ul}sd/"
python depmap_script.py -cf "${DATA_DIR}/19Q4/Achilles_gene_effect.csv" -cc "${OUTDIR}/_crispr_all_correlations.h5"\
  -rf "${DATA_DIR}/demeter/D2_combined_gene_dep_scores.csv" -rc "${OUTDIR}/_rnai_all_correlations.h5"\
  -b "${SIF_DIR}/belief_dict.pkl" --filter-type None --no-web-files\
  --explained-set /home/klas/repos/depmap_analysis/input_data/depmap/gene_sets/mito_function_depmap.csv gene\
  -lw "${SIF_DIR}/db_dump_df_lite.csv" -cstats 0.003076 0.056813 -rstats 0.006854 0.077614 -crange $ll $ul -rrange $ll $ul\
  -o "${OUTDIR}/${ll}_${ul}sd/" > "${OUTDIR}/${ll}_${ul}sd/stout.log"
echo "Finished running $ll-$ul SD"

ll=5;
echo "Running $ll-$ul SD"
echo "Output goes to ${OUTDIR}/${ll}_${ul}sd/"
python depmap_script.py -cf "${DATA_DIR}/19Q4/Achilles_gene_effect.csv" -cc "${OUTDIR}/_crispr_all_correlations.h5"\
  -rf "${DATA_DIR}/demeter/D2_combined_gene_dep_scores.csv" -rc "${OUTDIR}/_rnai_all_correlations.h5"\
  -b "${SIF_DIR}/belief_dict.pkl" --filter-type None --no-web-files\
  --explained-set /home/klas/repos/depmap_analysis/input_data/depmap/gene_sets/mito_function_depmap.csv gene\
  -lw "${SIF_DIR}/db_dump_df_lite.csv" -cstats 0.003076 0.056813 -rstats 0.006854 0.077614 -crange $ll $ul -rrange $ll $ul\
  -o "${OUTDIR}/${ll}_${ul}sd/" > "${OUTDIR}/${ll}_${ul}sd/stout.log"
echo "Finished running $ll-$ul SD"

