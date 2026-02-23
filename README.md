# arc_nocardia_omics
Code used to link genes involved in secondary metabolism to molecules of interest

# Files for further analysis:

### OneDrive data access link: [please email buultjensa@unimelb.edu.au to request access]
```
https://unimelbcloud-my.sharepoint.com/:f:/r/personal/buultjensa_unimelb_edu_au/Documents/arc_nocardia_omics?csf=1&web=1&e=gdEOui
```

### List of 128 isolates
```
fofn.txt
```

### Label vector files
```
/home/buultjensa/arc_nocardia_omics/ML/label_outputs/
```

### BiG-SCAPE feature table [n=1,053]
```
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/128_strain_by_family_matrix.tr.csv
```

### Panaroo feature table [n=150,747]
```
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98/gene_presence_absence_MATRIX.csv
```

### Panaroo expression checked feature table [n=150,747]
```
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98/128_panaroo_expression_checked_MATRIX.csv
```

### Report that link Panaroo and antiSMASH annotations with Panaroo and BiG-SCAPE clusters:
```
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/antismash_bgc_locus_report/128_antismash_bgc_locus_report_isolate_region_locus.csv

# Example:
$ head -2 128_antismash_bgc_locus_report_isolate_region_locus.csv
isolate,region_id,locus_tag,panaroo_gene_cluster,antismash_product_types,region_size_kb,bigscape_class,bigscape_family_number,bigscape_id
SP0005,SP0005_contig_1.region001,SP0005_00152,group_770,T1PKS,45.205,PKSI,4590,PKSI_4590

# Locus tag SP0005_00152 from isolate SP0005 is clustered into panaroo orthlog group group_770, is in antiSMASH region SP0005_contig_1.region001 of which is clustered into BiG-SCAPE gene cluster family PKSI_4590.
```

--------------------------------------------------------------------------

# Prokka command:
```
# Prokka output folders: /home/buultjensa/arc_nocardia_omics/prokka_gff/128

prokka --outdir ./Prokka_Output_240125/SP0129 \
--prefix SP0129 \
--locustag SP0129 \
/home/lsharkey/Projects/16_NocardiaGenomes/ARC_Strains/Jan_2025/Best_assemblies_compiled//SP0129.fasta \
--cpus 50 \
--genus nocardia
```

# Panaroo command:
```
# Panaroo output folders: /home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98

panaroo -i *.gff \
-o panaroo_test_c-0.98 \
--threads 16 \
--core_threshold 0.8 \
--clean-mode strict -c 0.98
```

# Check pangenome matrix for evidence of expression
```
python check_for_expression_v3.py \
<pangenome_csv> \
'<glob_pattern>' \
<output_csv>

# Command used for 301 media:
python check_for_expression_v3.py \
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98/gene_presence_absence.csv \
'/home/lsharkey/Projects/22_ARC_Project_Home/Proteomics/Full_dataset/outputs/output_080425/protein_tsvs/*301.tsv' \
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98/128_panaroo_expression_checked.csv
```

# antiSMASH command:
```
# Antismash output folders: /home/buultjensa/arc_nocardia_omics/copy_of_antismash/128

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

# BiG-SCAPE command:
```
# BiG-SCAPE output folders: /home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape

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
python make_strain_by_family_matrix_v3.py \
outfile.csv \
*/*clustering_c0.30.tsv

# Command used:
python make_strain_by_family_matrix_v3.py \
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/128_strain_by_family_matrix.csv \
/home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/*/*clustering_c0.30.tsv
```

# Make report that combines Prokka, panaroo, antiSMASH and BiG-SCAPE data
```
python antismash_bgc_locus_report.py \
gene_presence_absence.csv \
fofn.txt \
--antismash_dir [dir] \
--bigscape_run_dir [dir] \
--out_prefix [prefix]

# Command used:
python antismash_bgc_locus_report.py \
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/panaroo_test_c-0.98/gene_presence_absence.csv \
/home/buultjensa/arc_nocardia_omics/prokka_gff/128/128_fofn.txt \
--antismash_dir /home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/ \
--bigscape_run_dir /home/buultjensa/arc_nocardia_omics/copy_of_antismash/128/128_bigscape/network_files/2026-02-13_11-47-10_hybrids_glocal/ \
--out_prefix 128_antismash_bgc_locus_report
```











