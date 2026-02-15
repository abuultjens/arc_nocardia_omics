# arc_nocardia_omics
Code used to link genes involved in secondary metabolism to molecules of interest

# Files for further analysis:

### Labels
```
/home/buultjensa/arc_nocardia_omics/ML/label_outputs/
```

### BiG-SCAPE feature table
```
/home/buultjensa/arc_nocardia_omics/ML/128_strain_by_family_matrix.tr.csv
```

### Panaroo feature table
```
/home/buultjensa/arc_nocardia_omics/prokka_gff/154_gff/WO-SP0042/WO-subtree/WO-subtree_WO-647/panaroo_test_c-0.98/gene_presence_absence.csv
```

### Panaroo expression checked feature table
```
/home/buultjensa/arc_nocardia_omics/prokka_gff/154_gff/WO-SP0042/WO-subtree/WO-subtree_WO-647/panaroo_test_c-0.98/142_panaroo_expression_checked.csv
```

### Reports that link Panaroo and antiSMASH annotations with Panaroo and BiG-SCAPE clusters:
```
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/antismash_bgc_locus_report
```

--------------------------------------------------------------------------

# Prokka command:
```
prokka --outdir ./Prokka_Output_240125/SP0129 \
--prefix SP0129 \
--locustag SP0129 \
/home/lsharkey/Projects/16_NocardiaGenomes/ARC_Strains/Jan_2025/Best_assemblies_compiled//SP0129.fasta \
--cpus 50 \
--genus nocardia
```

# antiSMASH command:
```
antismash \
--genefinding-tool prodigal \
--cc-mibig \
--cb-general \
--cb-subcluster \
--cb-knownclusters \
--rre \
--cpus 6 \
/home/lsharkey/Projects/22_ARC_Project_Home/Genomics/prokka_annotated/Prokka_Output_240125/SP0585/SP0585.gbk \
--output-dir SP0585
```

# Panaroo command:
```
panaroo -i *.gff \
-o panaroo_test_c-0.98 \
--threads 16 \
--core_threshold 0.8 \
--clean-mode strict -c 0.98
```

# Check pangenome matrix for evidence of expression
```
/home/buultjensa/arc_nocardia_omics/prokka_gff/154_gff/WO-SP0042/WO-subtree/WO-subtree_WO-647/panaroo_test_c-0.98/check_for_expression_v2.py
check_for_expression_v2.py
```

# BiG-SCAPE command:
```
python3 /home/lsharkey/bigscape/BiG-SCAPE-1.1.9/bigscape.py \
-i /home/buultjensa/arc_nocardia_omics/copy_of_antismash \
-o bigscape-test-4 \
--pfam_dir /home/lsharkey/pfamdb \
--mibig \
--verbose \
-c 50 \
--clans-off
```

# Convert BiG-SCAPE output to a feature table:
```
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/make_strain_by_family_matrix_v2.py
```

# Make report that combines Prokka, panaroo, antiSMASH and BiG-SCAPE data
```
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/antismash_bgc_locus_report.py

python antismash_bgc_locus_report.py \ gene_presence_absence.csv \
fofn.txt \
--antismash_dir [dir] \
--bigscape_run_dir [dir] \
--out_prefix [prefix]
```











