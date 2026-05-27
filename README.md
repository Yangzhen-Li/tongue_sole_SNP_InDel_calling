You may need to adjust, according to the sequencing depth and sequencing technology. And can be easily adapted to be used in other species. 

cd /path/to/tongue_sole_variant_scripts

chmod +x *.sh

./01_fastp_qc.sh
./02_align_markdup_haplotypecaller.sh
./03_joint_genotyping_genomicsdb.sh
./04_select_snps_indels.sh
./05_hard_filter_snps_indels.sh
./06_vcftools_filter_wgs_snps.sh
./06b_vcftools_filter_20k_target_snps.sh
