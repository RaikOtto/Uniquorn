library("Uniquorn")

HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")

ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37" )

expect_that( ident_result, is_a("data.frame") )

expect_that( dim(ident_result)[1], equals( 60 )  )
expect_that( dim(ident_result)[2], equals( 9  ) )