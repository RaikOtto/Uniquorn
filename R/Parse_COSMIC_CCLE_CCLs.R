#' initiate_canonical_databases
#' 
#' Parses data into r list variable
#' 
#' @param cosmic_file The path to the Cosmic CLP file. The Cosmic file 
#' can be obtained from "https://cancer.sanger.ac.uk/cell_lines/download" and 
#' should be labeled "CosmicCLP_MutantExport.tsv.gz".
#' Ensure that the right reference genome is used
#' @param ccle_file The path to the ccle DNA genotype data file. 
#' It should be labeled "CCLE_mutations.csv".
#' Ensure that the right reference genome is used
#' @param ccle_sample_file The path to the CCLE sample file.
#' It should be labeled "sample_info.csv" containing both the
#' DepMap ID and corresponding cell line name.
#' @param ref_gen Reference genome version
#' @return Returns message if parsing process has succeeded
#' @import R.utils stringr
#' @usage
#' initiate_canonical_databases(
#'     cosmic_file = "CosmicMutantExport.tsv.gz",
#'     ccle_file = "CCLE_mutations.csv",
#'     ccle_sample_file = "sample_info.csv",
#'     ref_gen = "GRCH38"
#' )
#' @examples 
#' initiate_canonical_databases(
#'     cosmic_file = "CosmicMutantExport.tsv.gz",
#'     ccle_file = "CCLE_mutations.csv",
#'     ccle_sample_file = "sample_info.csv",
#'     ref_gen = "GRCH38"
#' )
#' @export
initiate_canonical_databases = function(
    cosmic_file = "CosmicMutantExport.tsv.gz",
    ccle_file = "CCLE_mutations.csv",
    ccle_sample_file = "sample_info.csv",
    ref_gen = "GRCH38"
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
    if (file.exists(ccle_file) && file.exists(ccle_sample_file)){
        message("Found CCLE and sample info: ", 
                file.exists(ccle_file) && file.exists(ccle_sample_file))
        parse_ccle_genotype_data(ccle_file, ccle_sample_file, ref_gen = ref_gen)
    }
    
    if ((!file.exists(cosmic_file)) & (!file.exists(ccle_file))){ 
        warning("Did neither find CCLE & CoSMIC file! Aborting.")
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
#' @importFrom stats aggregate
#' @importFrom data.table fread
#' @import IRanges
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
#' @keywords internal
parse_cosmic_genotype_data = function(cosmic_file, ref_gen = "GRCH38"){
    
    # Only read in columns specified with subset
    library_name = "COSMIC"
    package_path = system.file("", package = "Uniquorn")
    
    subset = c(5, 26)
    
    if (!grepl("CosmicMutantExport", cosmic_file)){ 
        stop("Warning. This is not the recommended COSMIC genotype file!",
            " The recommended file is the 'CosmicCLP_MutantExport.tsv.gz'",
            " file.")
        subset = c(5, 19)
    }
    
    cosmic_genotype_tab = data.table::fread(cosmic_file, select = subset,
       sep = "\t", showProgress = FALSE)
    colnames(cosmic_genotype_tab) = c("sample", "position") 
    
    # Remove CLs with NA or "" position
    empty_position <- c(which(cosmic_genotype_tab$position == ""), 
                        which(is.na(cosmic_genotype_tab$position)))
    cosmic_genotype_tab <- cosmic_genotype_tab[-empty_position,]
    
    # Extract and process coordinates and CL IDs
    message("Parsing Cosmic Coordinates")
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
    new_cls = data.table::data.table(
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
#' @param ccle_sample_file Path to CCLE sample file
#' @param ref_gen Reference genome version
#' @import IRanges
#' @importFrom stats aggregate
#' @importFrom data.table fread
#' @return The R Table sim_list which contains the CCLE fingerprints
#' @keywords internal
parse_ccle_genotype_data = function(ccle_file, ccle_sample_file, ref_gen = "GRCH38"){
    
    library_name = "CCLE"
    
    # Only read in columns specified with subset
    subset = c(4, 5, 6, 16)
    ccle_genotype_tab = data.table::fread(
        ccle_file,
        select = subset,
        sep = ",",
        showProgress = FALSE
    )
    
    subset = c(1,2)
    ccle_sample_tab = data.table::fread(
      ccle_sample_file,
      select = subset,
      sep = ",",
      showProgress = FALSE
    )
    
    depmap_ID_macthes <- match(ccle_genotype_tab$DepMap_ID, ccle_sample_tab$DepMap_ID)
    ccle_genotype_tab$Cell_line <- ccle_sample_tab$cell_line_name[depmap_ID_macthes]
    
    cls = ccle_genotype_tab[, gsub("\\_.*", "", Cell_line)] 
    cls = str_replace(cls, paste( "_", library_name, sep = "" ), "")
    
    coords = paste(
        ccle_genotype_tab$Chromosome,
        ccle_genotype_tab$Start_position,
        ccle_genotype_tab$End_position,
        sep = "_"
    )
    c_matches = match(coords, unique(coords), nomatch = 0)
    #new_cls <<- c()
    
    message("Aggregating CCLE CCL names")
    new_cls = data.table::data.table(
      "CLS" = cls,
      "Index" = c_matches
    )
    new_cls = new_cls[,lapply(.SD, paste, sep = "", collapse= ","), by = Index]
    
    # Extract and process coordinates and CL IDs
    g_mat = GenomicRanges::GRanges(
      seqnames = ccle_genotype_tab$Chromosome,
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
