context("Correct access and structure of stored database information")
names = read_library_names(ref_gen = "GRCH37")

test_that("show_contained_ccls returns correct object", {
    names = read_library_names(ref_gen = "GRCH37")
    ccls_all = show_contained_ccls()
  
    expect_is(ccls_all, "data.frame")
    expect_equal(dim(ccls_all)[2], 6)
    expect_equal(sort(unique(ccls_all$Library)), sort(names))
    expect_is(ccls_all[,1], "character")
    expect_is(ccls_all[,2], "character")

})

test_that("show_contained_variants_in_library returns correct object", {
    names = read_library_names(ref_gen = "GRCH37")
    if("CELLMINER" %in% names){
        
        g_ranges = show_contained_variants_in_library(
            library_name = "CELLMINER"
        )
        
        expect_is(g_ranges, "GRanges")
        expect_true(length(g_ranges) > 50000)
        expect_equal(sum(grepl("CELLMINER", g_ranges$Member_CCLs)), length(g_ranges))
    }
})

test_that("show_contained_variants_for_ccl returns the right object", {
    names = read_library_names(ref_gen = "GRCH37")
    ccls_all = show_contained_ccls()
    if("CELLMINER" %in% names){
      
        index = ccls_all$Library == "CELLMINER"
        ccls_all = ccls_all[index,]
        
        ccl_index = sample(1:length(ccls_all), 1, replace = FALSE)
        ccl = ccls_all$CCL[ccl_index]
        
        variants = show_contained_variants_for_ccl(
            name_ccl = ccl, library_name = "CELLMINER"
        )
  
        expect_is(variants, "GRanges")
        expect_true(length(variants) > 1)
        
        ccl_index = sample(1:length(ccls_all), 1, replace = FALSE)
        ccl = ccls_all$CCL[ccl_index]
        
        variants = show_contained_variants_for_ccl(
          name_ccl = ccl, library_name = "CELLMINER"
        )
        
        expect_is(variants, "GRanges")
        expect_true(length(variants) > 1)
  }
})