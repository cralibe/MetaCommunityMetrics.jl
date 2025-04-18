function average_MVNH_dissimilarity_v2(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; var_names::Vector{String}=String[])
    @assert length(presence_absence) == nrow(data) "Mismatch in data length"
    @assert length(species) == nrow(data) "Mismatch in species length"

    # Pre-compute unique species list once
    unique_species = unique(species)
    
    # Convert DataFrame to matrix once (this is much more efficient)
    data_matrix = Matrix(data)
    
    # Pre-calculate indices for each species with presence > 0
    species_indices = Dict{String, Vector{Int}}()
    for (i, (sp, pres)) in enumerate(zip(species, presence_absence))
        if pres > 0
            if !haskey(species_indices, sp)
                species_indices[sp] = Int[]
            end
            push!(species_indices[sp], i)
        end
    end
    
    # Initialize final result
    total_dissimilarity = 0.0
    pair_count = 0
    
    # Process each species pair exactly once
    for i in 1:length(unique_species)-1
        sp1 = unique_species[i]
        
        # Skip if no indices for this species
        if !haskey(species_indices, sp1) || isempty(species_indices[sp1])
            continue
        end
        
        # Get data for species 1 once
        indices1 = species_indices[sp1]
        if isempty(indices1)
            continue
        end
        
        # Convert to DataFrame once for this species
        filtered_df_1 = DataFrame(view(data_matrix, indices1, :), var_names)
        
        for j in (i+1):length(unique_species)
            sp2 = unique_species[j]
            
            # Skip if no indices for this species
            if !haskey(species_indices, sp2) || isempty(species_indices[sp2])
                continue
            end
            
            # Get data for species 2
            indices2 = species_indices[sp2]
            if isempty(indices2)
                continue
            end
            
            # Convert to DataFrame for species 2
            filtered_df_2 = DataFrame(view(data_matrix, indices2, :), var_names)
            
            # Compute the MVNH_dissimilarity result
            result = MVNH_dissimilarity(filtered_df_1, filtered_df_2; var_names)
            
            # Filter the result DataFrame to get the Bhattacharyya distance
            filtered_result = filter(row -> row[:metric] == "Bhattacharyya_distance", result)
            bhattacharyya_distance = filtered_result[1, :total]
            
            # Add to total
            total_dissimilarity += bhattacharyya_distance
            pair_count += 1
        end
    end
    
    # Return average
    return pair_count > 0 ? total_dissimilarity / pair_count : 0.0
end