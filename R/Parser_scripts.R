parse_cosmic_genotype_data = function( cosmic_file ){
    
    exclude_cols = c(rep("NULL",4),"character",rep("NULL",12),"character","character",rep("NULL",13))
    
    #cosmic_genotype_file = "/Users/raik_000/Dropbox/PhD/Uniquorn_project//Raw_data/CosmicCLP_MutantExport.tsv"
    cosmic_genotype_tab = read.table(cosmic_genotype_file, sep ="\t", header =T, fill = T, col.names = exclude_cols)
    
    col_names_sim_list =  
    sim_list = matrix()
    
}

parse_ccle_genotype_data = function( ccle_file ){
    
    exclude_cols = c(rep("NULL",4),"character",rep("NULL",12),"character","character",rep("NULL",13))
    
    #cosmic_genotype_file = "/Users/raik_000/Dropbox/PhD/Uniquorn_project//Raw_data/CosmicCLP_MutantExport.tsv"
    cosmic_genotype_tab = read.table(cosmic_genotype_file, sep ="\t", header =T, fill = T, col.names = exclude_cols)
    
    col_names_sim_list =  
        sim_list = matrix()
    
}

