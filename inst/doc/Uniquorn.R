## ------------------------------------------------------------------------
library("Uniquorn")
HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37"  )
head( ident_result )

