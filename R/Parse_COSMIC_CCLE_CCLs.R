#' initiate_canonical_databases
#' 
#' Parses data into r list variable
#' 
#' @param cosmic_file The path to the cosmic DNA genotype data file. 
#' Ensure that the right reference genome is used
#' @param ccle_file The path to the ccle DNA genotype data file. 
#' Ensure that the right reference genome is used
#' @param ref_gen Reference genome version
#' @return Returns message if parsing process has succeeded
#' @import R.utils stringr
#' @usage 
#' initiate_canonical_databases(
#'     cosmic_file = "CosmicCLP_MutantExport.tsv",
#'     ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
#'     ref_gen = "GRCH37"
#' )
#' @examples 
#' initiate_canonical_databases(
#'     cosmic_file = "CosmicCLP_MutantExport.tsv",
#'     ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
#'     ref_gen = "GRCH37"
#' )
#' @export
initiate_canonical_databases = function(
    cosmic_file = "CosmicCLP_MutantExport.tsv",
    ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
    ref_gen = "GRCH37"
){

    message("Reference genome: ", ref_gen)

    # Parse CoSMIC file
    if (file.exists(cosmic_file)){
        message("Found CoSMIC: ", file.exists(cosmic_file))
        
        if (grepl(".gz$", cosmic_file, ignore.case = TRUE)){
            gunzip(cosmic_file, overwrite = TRUE)
            cosmic_file = gsub(".gz$", "", cosmic_file, ignore.case = TRUE)
        }
        
        parse_cosmic_genotype_data( cosmic_file, ref_gen = ref_gen )
    }

    # Parse CCLE file
    if (file.exists(ccle_file)){
        message("Found CCLE: ", file.exists(ccle_file))
        parse_ccle_genotype_data(ccle_file, ref_gen = ref_gen)
    }
    
    if ((!file.exists(cosmic_file)) & (!file.exists(ccle_file))){ 
        warning("Did neither find CCLE & CoSMIC CLP file! Aborting.")
    } else {
        message("Finished parsing, ",
            "aggregating over parsed Cancer Cell Line data")
        
        message("Finished aggregating, saving to database")

        message("Initialization of Uniquorn DB finished")
    }
}

#' parse_cosmic_genotype_data
#' 
#' Parses cosmic genotype data
#' 
#' @param cosmic_file Path to cosmic clp file in hard disk
#' @param ref_gen Reference genome version
#' @import data.table
#' @importFrom IRanges IRanges
#' @importFrom stats aggregate
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
parse_cosmic_genotype_data = function(cosmic_file, ref_gen = "GRCH37"){
    
    # Only read in columns specified with subset
    library_name = "COSMIC"
    package_path = system.file("", package = "Uniquorn")
    
    subset = c(5, 24)
    
    if (!grepl("CosmicCLP_MutantExport", cosmic_file)){ # MutantExport
        stop("Warning. This is not the recommended COSMIC genotype file!",
            " The recommended file is the 'CosmicCLP_MutantExport.tsv.gz'",
            " file.")
        subset = c(5, 19)
    }
    
    cosmic_genotype_tab = fread(cosmic_file, select = subset,
        sep = "\t", showProgress = FALSE)
    colnames(cosmic_genotype_tab) = c("sample", "position")
    
    # Extract and process coordinates and CL IDs
    message("Parsing Cosmic Coordinates, that might take some time")
    coords = cosmic_genotype_tab[, gsub(":|-", "_", position)]
    seq_name = vapply(strsplit(coords, "_"), `[`, 1, FUN.VALUE = character(1))
    starts = vapply(strsplit(coords, "_"), `[`, 2, FUN.VALUE = character(1))
    ends = vapply(strsplit(coords, "_"), `[`, 3, FUN.VALUE = character(1))

    cls = cosmic_genotype_tab[
        ,gsub("/|(|])| ", "", sample, ignore.case = TRUE)]
    cls[cls == "KM-H2"] = "KMH2"
    cls[cls == "KMH-2"] = "KMH2ALTERNATIVE"
    
    c_matches = match(coords, unique(coords), nomatch = 0)
    
    message("Aggregating Cosmic CCL names")
    new_cls = data.table(
      "CLS" = cls,
      "Index" = c_matches
    )
    new_cls = new_cls[,lapply(.SD, paste, sep = "", collapse= ","), by = Index]
    
    # Extract and process coordinates and CL IDs
    g_mat = GenomicRanges::GRanges(
        seqnames = seq_name,
        IRanges(
            start = as.integer( starts ),
            end = as.integer( ends )
        )
    )
    g_mat = unique(g_mat)
    mcols(g_mat)$Member_CCLs = new_cls$CLS
    mcols(g_mat)$Member_CCLs = str_replace_all(mcols(g_mat)$Member_CCLs, 
        pattern = paste( "_", library_name, sep =""), "")
    
    message("Normalizing CCL names")
    
    write_mutation_grange_objects(
        g_mat = g_mat,
        library_name = library_name,
        ref_gen = ref_gen, 
        mutational_weight_inclusion_threshold = 0
    )
    
    write_w0_and_split_w0_into_lower_weights(
        g_mat = g_mat,
        ref_gen = ref_gen,
        library_name = "COSMIC"
    )
    message("Finished parsing Cosmic")
}

#' parse_ccle_genotype_data
#' 
#' Parses ccle genotype data
#' 
#' @param ccle_file Path to CCLE file on hard disk
#' @param ref_gen Reference genome version
#' @import data.table
#' @importFrom IRanges IRanges
#' @importFrom stats aggregate
#' @return The R Table sim_list which contains the CCLE fingerprints
parse_ccle_genotype_data = function(ccle_file, ref_gen = "GRCH37"){
    
    library_name = "CCLE"
    
    # Only read in columns specified with subset
    subset = c(5, 6, 7, 16)
    ccle_genotype_tab = data.table::fread(
        ccle_file,
        select = subset,
        sep = "\t",
        showProgress = FALSE
    )
    
    cls = ccle_genotype_tab[, gsub("\\_.*", "", Tumor_Sample_Barcode)]
    cls = str_replace(cls, paste( "_", library_name, sep = "" ), "")
    
    coords = paste(
        str_replace(ccle_genotype_tab$Chromosome, pattern = "^chr", "" ),
        ccle_genotype_tab$Start_position,
        ccle_genotype_tab$End_position,
        sep = "_"
    )
    c_matches = match(coords, unique(coords), nomatch = 0)
    #new_cls <<- c()
    
    message("Aggregating CCLE CCL names")
    new_cls = data.table(
      "CLS" = cls,
      "Index" = c_matches
    )
    new_cls = new_cls[,lapply(.SD, paste, sep = "", collapse= ","), by = Index]
    
    # Extract and process coordinates and CL IDs
    g_mat = GenomicRanges::GRanges(
        seqnames = str_replace(
            ccle_genotype_tab$Chromosome, pattern = "^chr", "" ),
        IRanges(
            start = ccle_genotype_tab$Start_position,
            end = ccle_genotype_tab$End_position
        )
    )
    g_mat = unique(g_mat)
    mcols(g_mat)$Member_CCLs = new_cls$CLS
    
    # write the stats part
    write_w0_and_split_w0_into_lower_weights(
        g_mat = g_mat,
        ref_gen = ref_gen,
        library_name = "CCLE"
    )
    
    message("Finished parsing CCLE")
}