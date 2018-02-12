context("Correct add and remove of custom Cancer Cell Line")

test_that("adding and removing a custom Cancer Cell Line works succesfully", {

    # Check add
    HT29_vcf_file = system.file("extdata/HT29_TEST.vcf", package = "Uniquorn")
    expect_message(
        add_custom_vcf_to_database(
            vcf_input_files = HT29_vcf_file,
            ref_gen = "GRCH37",
            library_name = "CELLMINER",
            n_threads = 1,
            test_mode = TRUE
        ),
        "Finished"
    )
  
    #ccls_all = show_contained_ccls()
    #variants = show_contained_variants_for_ccl(name_ccl = "HT29_TEST", library_name = "CELLMINER")
  
    #expect_equal(dim(ccls_all)[1], ccls_before + 1)
    #expect_true(any(grepl("HT29_TEST", ccls_all[,1], fixed = TRUE)))
    #expect_equal(length(variants), 24969)
  
    # Check remove
    expect_message(
        remove_ccls_from_database(
            ccl_names = "786_0",
            ref_gen = "GRCH37",
            library_name = "CELLMINER",
            test_mode = TRUE
        ),
        "Finished removing all cancer cell lines"
    )
    #ccls_all = show_contained_ccls()
    #variants = show_contained_variants_for_ccl(name_ccl = "HT29_TEST", library_name = "CELLMINER")
  
    #expect_equal(dim(ccls_all)[1], ccls_before)
    #expect_false(any(grepl("HT29_TEST", ccls_all[,1], fixed = TRUE)))
    #expect_equal(length(variants), 0)
})