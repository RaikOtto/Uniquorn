context("Correct identification results for test data")

test_that("identification of test data returns correct results", {
    HT29_vcf_file = system.file("extdata", "HT29.vcf.gz", package = "Uniquorn")
    ident_result = identify_vcf_file(HT29_vcf_file, ref_gen = "GRCH37", verbose = FALSE)
    
    libraries = read_library_names(ref_gen = "GRCH37")
    ccls_all = show_contained_ccls()
    ccls = ccls_all[,1]
    ccls = str_replace_all(ccls, pattern = 
        paste(c("\\_","\\-","\\(","\\)"), collapse = "|", sep = ""), "")
    ccls = str_to_upper(ccls)
    occ = sum(ccls %in% "HT29")
    index = which(ccls %in% "HT29")
    names = ccls_all[index, 1]
    
    expect_is(ident_result, "data.frame")
    expect_equal(dim(ident_result)[2], 7)
  
    index = which(ident_result[,1] %in% names)
    expect_equal(length(index), occ)
    expect_true(all(ident_result$P_value_sig[index]))
    expect_true(all(ident_result$Above_Penality[index]))
    expect_true(all(ident_result$Identification_sig[index]))
    expect_equal(sort(unique(as.character(ident_result$Library))), sort(libraries))

})