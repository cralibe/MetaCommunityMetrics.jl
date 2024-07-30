# src/Internal.jl

module Internal

using DataFrames

# Internal function
# These functions is for internal use and supports the public functions.

## Internal fuctions for DNCI.jl
#A function assign single site to the nearest group
function assign_single_site_to_nearest_group(grouped_df, single_site, min_group)
    # Extract latitude and longitude of the single site
    single_lat = single_site.Latitude[1]
    single_lon = single_site.Longitude[1]
    
    min_distance = Inf
    nearest_group_id = nothing

    # Iterate over other groups
    for row in eachrow(grouped_df[grouped_df.Group .!= min_group, :])
        dist = sqrt((single_lat - row.Latitude)^2 + (single_lon - row.Longitude)^2)
        if dist < min_distance
            min_distance = dist
            nearest_group_id = row.Group
        end
    end

    if isnothing(nearest_group_id)
        error("No valid nearest group found for the single site.")
    end

    # Assign the single site to the nearest group
    #single_site.Group[1] = nearest_group_id

    return nearest_group_id
end

#A function that find the nearest site to the group with the fewest sites
function find_nearest_site(grouped_df, min_group)
    # Extract coordinates of sites in the group with the fewest species
    sites_in_min_group = grouped_df[grouped_df.Group .== min_group, :]
    
    # Check if there are any sites in the minimum group
    if isempty(sites_in_min_group)
        error("No sites found in the minimum group.")
    end

    # Calculate centroid of the cluster
    centroid_latitude = mean(sites_in_min_group.Latitude)
    centroid_longitude = mean(sites_in_min_group.Longitude)
    
    # Initialize variables to keep track of the nearest site and distance
    min_distance = Inf
    site_with_minimum_distance = nothing  # Initialize with nothing to handle cases where no site is found

    # Filter for other groups
    other_group = grouped_df[grouped_df.Group .!= min_group, :]

    # Check if there are other groups available
    if isempty(other_group)
        error("No other groups found to compare.")
    end

    # Iterate over each site in other_group
    for row in eachrow(other_group)
        # Calculate the Euclidean distance between the centroid and the current site
        dist = sqrt((centroid_latitude - row.Latitude)^2 + (centroid_longitude - row.Longitude)^2)
     
        # Update the nearest site if this site is closer
        if dist < min_distance
            min_distance = dist
            site_with_minimum_distance = row
        end
    end

    # Check if a nearest site was found
    if isnothing(site_with_minimum_distance)
        error("No nearest site could be determined.")
    end
        
    return site_with_minimum_distance
end
#A function that checks the condition for the groupings only
function check_conditions(subset_df::DataFrame)
    unique_groups = unique(subset_df.Group)
    sites_per_group = Dict(g => sum(subset_df.Group .== g) for g in unique_groups)
    richness_per_group = Dict(g => sum(subset_df.Total_Richness[subset_df.Group .== g]) for g in unique_groups)
    for i in unique_groups
        for j in unique_groups
            if i != j
                site_diff_ratio = abs(sites_per_group[i] - sites_per_group[j]) / max(sites_per_group[i], sites_per_group[j])
                species_diff_ratio = abs(richness_per_group[i] - richness_per_group[j]) / max(richness_per_group[i], richness_per_group[j])
                if site_diff_ratio > 0.3 || species_diff_ratio > 0.4
                    return false
                end
            end
        end
    end

    if any(value < 5 for value in values(sites_per_group)) || length(unique_groups) <= 1
        return false
    end

    return true
end

#A function that checks the condition for the groupings and fix the groupings
function check_condition_and_fix(grouped_df)
    max_iterations = 1000  # Define a maximum number of iterations to prevent infinite loops
    iteration = 0

    while iteration < max_iterations
        unique_groups = unique(grouped_df.Group)
        sites_per_group = Dict(g => sum(grouped_df.Group .== g) for g in unique_groups)
        richness_per_group = Dict(g => sum(grouped_df.Total_Richness[grouped_df.Group .== g]) for g in unique_groups)

        condition_met = true

        # Checking conditions including group size constraints within main condition checks
        for i in unique_groups
            for j in unique_groups
                if i != j
                    site_diff_ratio = abs(sites_per_group[i] - sites_per_group[j]) / max(sites_per_group[i], sites_per_group[j])
                    species_diff_ratio = abs(richness_per_group[i] - richness_per_group[j]) / max(richness_per_group[i], richness_per_group[j])                
                
                    if site_diff_ratio > 0.3 || species_diff_ratio > 0.4 
                        condition_met = false
                        break
                    end
                end
                
                if !condition_met
                    break
                end
            end
            
            if !condition_met
                break
            end
        end

        if any(value < 5 for value in values(sites_per_group)) || length(unique_groups) <= 1
            condition_met = false
        end
        
        if condition_met
            println("Condition met, stopping iteration.")
            break
        else
             
            #println("Condition not met, adjusting data...") 

            min_group = first(argmin(sites_per_group))

            # Check if there is only one site in the minimum group
            sites_in_min_group = grouped_df[grouped_df.Group .== min_group, :]
            if nrow(sites_in_min_group) == 1
            # Find the nearest group to this single site
            grouped_df.Group[grouped_df.Group .== min_group, :] .= assign_single_site_to_nearest_group(grouped_df, sites_in_min_group, min_group)
            else
            site_with_minimum_distance = find_nearest_site(grouped_df, min_group) #Find the nearest site to the group with the fewest sites
            grouped_df.Group[grouped_df.Patch .== site_with_minimum_distance.Patch] .= min_group #assign the nearest site to the group with the fewest sites
            end
        end

        iteration += 1
    end

    if iteration == max_iterations
        println("Reached maximum iterations without meeting condition.")
    end

    return grouped_df
end

end