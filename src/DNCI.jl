# src/DNCI.jl

using ..Internal


"""
    create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}, total_richness::Vector{Int}) -> Dict{Int, DataFrame}

Create clusters for each unique time step in a dataset. Only presnece-absence data can be used.

# Arguments
- `time::Vector`: A vector indicating the time each sample was taken.
- `latitude::Vector`: A vector indicating the latitude of each sample.
- `longitude::Vector`: A vector indicating the longitude of each sample.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.
- `total_richness::Vector`: A vector indicating the total species richness at each plot at each time step.

# Returns
- `Dict{Int, DataFrame}`: A dictionary where each key is a unique time from the dataset and each value is a DataFrame for that time with an added `Group` column indicating the assigned cluster.

# Details
This function performs hierarchical clustering on geographical coordinates for each unique time step. It aims to balance the clusters based on the number of sites and species richness, ensuring that no group has less than five sites and that there are at least two groups. If conditions for clustering balance are not met (like groups having less than five sites or only one group), it iteratively adjusts the clusters by reassigning sites to improve group balance.

"""
function create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}, total_richness::Vector{Int})
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Patch = patch, Total_Richness = total_richness)
    
    for t in unique(df[:, :Time])
        subset_df = filter(row -> row[:Time] == t, df)
        coordinates = select(subset_df, [:Latitude, :Longitude])
        distances = Distances.pairwise(Euclidean(), Matrix(coordinates), dims=1)
        
        num_clusters = max(div(length(unique(subset_df.Patch)), 5), 2) # Ensure that the number of clusters is at least 2 and each group has at least 5 sites
        condition_met = false

        while !condition_met && num_clusters > 1
            agglo_result = hclust(distances, linkage=:complete)
            assignments = cutree(agglo_result, k=num_clusters)
            subset_df.Group = assignments
            subset_df = Internal.check_condition_and_fix(subset_df)
            condition_met = Internal.check_conditions(subset_df)
            
            if !condition_met
                num_clusters -= 1
                if num_clusters < 2
                    println("Warning: Cluster count fell below 2, which is not permissible for clustering. Groups assigned as missing.")
                    subset_df.Group .= missing
                end
            end
        end

        grouping_dict[t] = subset_df
    end
    if length(grouping_dict) != length(unique(df.Time))
        println("Warning: Some time steps are missing!!!")
    end    
    return grouping_dict
end

"""
    plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::Union{AbstractVector, String}) 

Plots clusters of data points on a scatter plot using their geographic coordinates and cluster assignments.

# Arguments
- `latitude::Vector{Float64}`: A vector of latitude coordinates.
- `longitude::Vector{Float64}`: A vector of longitude coordinates.
- `group::Union{AbstractVector, String}`: A vector or string indicating the cluster assignments for each data point.

# Details
- The function assigns a unique color to each cluster and plots the points based on their geographic coordinates.
- The points are colored according to their cluster assignment.
- The plot includes black borders around the markers for better visibility.

# Example
```julia
latitudes = [35.0, 35.0, 35.5, 35.5, 35.5, 36.0, 36.0]
longitudes = [-110.0, -109.5, -109.5, -109.0, -108.0, -109.5, -108.0]
groups = [1, 1, 1, 1, 2, 1, 2]

plot_clusters(latitudes, longitudes, groups)
```
"""
function plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::Union{AbstractVector, String})
    # Get unique cluster IDs and assign numeric identifiers
    unique_clusters = unique(group)
    cluster_map = Dict(cluster => i for (i, cluster) in enumerate(unique_clusters))

    # Assign numeric identifiers to cluster IDs
    numeric_ids = [cluster_map[cluster] for cluster in group]

    # Define color palette for clusters
    colors = distinguishable_colors(length(unique_clusters), colorant"blue")

    # Plot the points, color by cluster ID
    scatter(longitude, latitude, marker_z=numeric_ids,
            xlabel="Longitude", ylabel="Latitude", title="Cluster Visualization",
            legend=false, markerstrokecolor=:black, markerstrokewidth=0.5,
            markersize=5, color=colors[numeric_ids], label=false)
end


"""
    DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000, count::Bool=true) -> DataFrame

Calculates the dispersal-niche continuum index (DNCI) for multiple groups, a metric proposed by Vilmi(2021) (doi: 10.1111/ecog.05356).

# Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites or samples) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `count::Bool=true`: A flag indicating whether the numeber of permutations is printed. Default is `false`.

# Returns
- `DataFrame`: A DataFrame containing the DNCI results for each pair of groups.


# Examples
```julia
# Example usage of DNCI_multigroup
comm = [1 0 0 1 0;
        1 1 0 0 0;
        0 1 1 0 0;
        0 0 1 1 1;
        1 0 0 0 1;
        0 1 1 0 1]

groups = ["A", "A", "B", "B", "C", "C"]
Nperm = 1000
count = true

result = DNCI_multigroup(comm, groups, Nperm, count)
println(result)
```
"""
function DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; count::Bool=true) #for presence-absence data only
    group_combinations = collect(combinations(unique(sort(groups)),2))

    ddelta = DataFrame()

    for i in 1:size(group_combinations,1)
        # Create an empty dictionary to hold the split data
        splitx = Dict()

        # Assume comm is a matrix and groups is a vector indicating the group for each row in comm
        for row in 1:size(comm, 1)  # Iterate over the rows of the matrix
            group = groups[row]  # Identify the group of the current row
            current_row = comm[row, :]   # Extract the entire row as a vector, make it a 1-row matrix
            row_matrix = reshape(current_row, 1, length(current_row))  # Reshape row to be a 1-row matrix

            # Check if the group already exists in the dictionary
            if haskey(splitx, group)
            # Vertically concatenate the new row matrix with the existing matrix for this group
            splitx[group] = vcat(splitx[group], row_matrix)
            else
            # Initialize the group with the current row as a 1-row matrix
                splitx[group] = row_matrix
            end
        end

        # Safely access and concatenate matrices from two groups
        if haskey(splitx, group_combinations[i][1]) && haskey(splitx, group_combinations[i][2])
            paired_x = vcat(splitx[group_combinations[i][1]], splitx[group_combinations[i][2]])
        else
            error("One of the groups does not exist in the dictionary.")
        end
        # Calculate a logical array indicating non-zero sum columns
        non_zero_sum_columns = sum(paired_x, dims=1) .!= 0
        # Convert logical index to actual column indices
        column_indices = findall(x -> x, non_zero_sum_columns[:]) 
    
        # Subset the matrix using these indices
        paired_x = paired_x[:, column_indices]

        group_pair = vcat(
        fill(group_combinations[i][1], size(splitx[group_combinations[i][1]], 1)),
        fill(group_combinations[i][2], size(splitx[group_combinations[i][2]], 1)))

        DNCI_result = Internal.DNCI_ses(paired_x, group_pair, Nperm; count)

        append!(ddelta, DNCI_result, promote=true)
    end
    return ddelta
end
