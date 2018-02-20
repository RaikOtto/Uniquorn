#' add_p_q_values_statistics
#' 
#' A hypergeometric distribution-assumption allows to calculate the 
#' p-values for a significant or non-significant overlap in this function
#' 
#' \code{add_p_q_values_statistics} Calculates the p-values
#' 
#' @param g_query IRanges object that contains the query variants
#' @param match_t A table that contains the nubmber of matching variants
#' @param p_value Threshold for the significance p-value
#' @param ref_gen Reference genome version
#' @param minimum_matching_mutations Manual lower amount of matching
#' mutations require for a significant match between a query and a 
#' reference
#' @param top_hits_per_library limits significant similarities to the first 
#' n hits
#' @import stringr
#' @importFrom stats pbinom
#' @return R table with a statistic 
add_p_q_values_statistics = function(
    g_query,
    match_t,
    p_value,
    ref_gen,
    minimum_matching_mutations,
    top_hits_per_library
){
    
    library_names = read_library_names(ref_gen = ref_gen)
    p_values = rep(1.0, nrow(match_t))
    balls_in_query = length(g_query$Member_CCLs)
    sig_vec <<- c()
    
    for (library_name in library_names){
        
        index_library = which( match_t$Library == library_name  )
        white_balls_possible = as.integer(as.character(
            match_t$All_variants))[index_library]
        
        if( length(white_balls_possible) == 0){
            stop("There is only one CL in the custom set.",  
                "A confidence score calculation is only ", 
                "possible with more than one sample!")
        }
        
        white_balls_found = as.integer(
            as.character(match_t$Matches))[index_library]
        white_balls_found_least_one = sum(white_balls_found > 0)
        black_balls = sum(white_balls_possible) - white_balls_possible
        
        background_cls_traces = sum(
            white_balls_found >= mean(white_balls_found)
        )
        
        #likelihood_found = ( white_balls_possible / 
        #    sum(white_balls_possible) ) **
        #    (white_balls_found / sum(white_balls_found))
        likelihood_found = white_balls_possible / sum(white_balls_possible)

        q = white_balls_found - 1
        q[q < 0]   = 0
        
        p_values_panel = as.double(stats::pbinom(
            q = q,
            size = sum(white_balls_found),
            #size = balls_in_query,
            p = likelihood_found,
            lower.tail = FALSE 
        ))
        
        p_values_panel[white_balls_found == white_balls_possible] = 0.0
        p_values_panel[p_values_panel >= 1.0] = 1.0
        p_values_panel[white_balls_found == 0] = 1.0
        p_values_panel[white_balls_found < minimum_matching_mutations] = 1.0
        
        p_values[index_library] = p_values_panel
        n_hits_vec = order( as.double( p_values[index_library] ) )
        n_hits_vec[as.integer(n_hits_vec) >
            as.integer(top_hits_per_library)] = FALSE
        n_hits_vec[as.integer(n_hits_vec) > 0 ] = TRUE
        n_hits_vec = as.logical(n_hits_vec)
        
        sig_vec <<- c(sig_vec, n_hits_vec)
    }
    
    match_t$P_values = p_values
    match_t$P_value_sig = match_t$P_values <= p_value
    match_t$P_value_sig = match_t$P_value_sig & sig_vec

    return(match_t)
}

#' add_penalty_statistics
#' 
#' Add penalty statistics to results
#' 
#' @param match_t object that contains the matching variants
#' @param minimum_matching_mutations a numerical giving the minimum amount of
#'  mutations that has to match between query and training sample for a 
#'  positive prediction
#' @importFrom stats integrate
#' @importFrom stats pbeta
#' @return The updated statistics
add_penality_statistics = function(match_t, minimum_matching_mutations){

    matching_variants = as.integer(match_t$Matches)
    minimum_matching_mutations = min(matching_variants)
    matching_variants = matching_variants[ matching_variants > 0]
    
    if ( length(matching_variants) == 0 ){
        
        penalty = 0
        penalty_mutations = 0
        
    } else{ 
        
        mean_match = mean( matching_variants )
        max_match  = max( matching_variants )
        
        penalty = integrate(
            f = pbeta,
            0,
            1,
            max_match,
            max_match/ mean_match,
            stop.on.error = FALSE
            )$value
            
        penalty_mutations = ceiling(
            mean_match + (max_match * penalty) / (1 - penalty)
        )
    }
    
    if ( ( minimum_matching_mutations == 0 ) & ( penalty != 0.0) ){
        
        message(
            "Correcting the background due to traces of random, ", 
            "scale-freeness amounts of matches, requiring at least ", 
            as.character(penalty_mutations), " variants to match."
        )
    }
    match_t$Above_Penality = 
        as.integer(as.character(match_t$Matches)) > penalty
    
    return(match_t)
}