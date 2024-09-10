# src/DNCI.jl

using ..Internal


"""
    create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}, total_richness::Vector{Int}) -> Dict{Int, DataFrame}

This function creates clusters (groupings of patches/sites) for each unique time step in a dataset which can then used for calculating DNCI. Only presnece-absence data can be used. Please remove singletons (taxa/species that occuring at one patch/site within a time step) before using this function.

Arguments
- `time::Vector`: A vector indicating the time each sample was taken.
- `latitude::Vector`: A vector indicating the latitude of each sample.
- `longitude::Vector`: A vector indicating the longitude of each sample.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample. At least 10 patches are required for clustering.
- `total_richness::Vector`: A vector indicating the total species richness at each plot at each time step.

Returns
- `Dict{Int, DataFrame}`: A dictionary where each key represents a unique time point from the input data, with the corresponding value being a `DataFrame` for that time step. Each `DataFrame` contains the following columns: `Time`, `Latitude`, `Longitude`, `Patch`, `Total_Richness`, and `Group` (indicating the assigned cluster).

Details
This function performs hierarchical clustering on the geographical coordinates of sampling patches/sites at each unique time step, assuming that organism dispersal occurs within the study region. It incorporates checks and adjustments to ensure the following conditions are met: at least 2 clusters, a minimum of 5 patches/sites per cluster, and that the variation in the number of taxa/species and patches/sites per group does not exceed 40% and 30%, respectively. These conditions are critical for calculating an unbiased DNCI value, and the function will issue warnings if any are not fulfilled.

Example
```jildoctest
julia> using DataFrames, CSV, MetaCommunityMetrics, Random, Distributions

julia> path = joinpath(pkgdir(MetaCommunityMetrics), "data

julia> time = vcat(fill(1, 20), fill(2, 20), fill(3, 20)) 
60-element Vector{Int64}:
 1
 1
 1
 ⋮
 3
 3

julia> latitude = rand(Uniform(34.0, 35.0), 60)
60-element Vector{Float64}:
 34.57986212013413
 34.41129411794985
 34.97213608245547
  ⋮
 34.57432348527832
 34.67764990759958

julia> longitude = rand(Uniform(-118.0, -117.0), 60)
60-element Vector{Float64}:
 -117.91655991056788
 -117.47420433613088
 -117.15935908052177
    ⋮
 -117.92682908601348
 -117.49824402144651

julia> patch = rand(1:10, 60) 
60-element Vector{Int64}:
 10
  2
  4
  ⋮
  6
  9

julia> total_richness = rand(5:50, 60)
60-element Vector{Int64}:
 37
 40
  6
  ⋮
 48
 28

julia> result = create_clusters(time, latitude, longitude, patch, total_richness)
Dict{Int64, DataFrame} with 3 entries:
  2 => 20×6 DataFrame…
  3 => 20×6 DataFrame…
  1 => 20×6 DataFrame…

julia> println(result[1])
20×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1   34.5799   -117.917     10              37      1
   2 │     1   34.4113   -117.474      2              40      2
   3 │     1   34.9721   -117.159      4               6      2
   4 │     1   34.0149   -117.476      9              45      1
   5 │     1   34.5204   -117.987      7              29      1
   6 │     1   34.6396   -117.594      5              18      2
   7 │     1   34.8396   -117.876     10               7      1
   8 │     1   34.9671   -117.351      1              34      2
   9 │     1   34.7898   -117.085      4              23      2
  10 │     1   34.696    -117.666      9              24      1
  11 │     1   34.5667   -117.652      8              19      1
  12 │     1   34.5364   -117.471     10              45      2
  13 │     1   34.7114   -117.042      5              43      2
  14 │     1   34.1039   -117.979      2              46      1
  15 │     1   34.8067   -117.637     10               8      2
  16 │     1   34.8705   -117.997      7              33      1
  17 │     1   34.9627   -117.369      9              48      1
  18 │     1   34.1512   -117.165      4              27      2
  19 │     1   34.7154   -117.12       8              34      1
  20 │     1   34.9395   -117.247      5              13      2

```
"""
function create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}, total_richness::Vector{Int})
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Patch = patch, Total_Richness = total_richness)
    
    for t in unique(df[:, :Time])
        subset_df = filter(row -> row[:Time] == t, df)
        num_sites = nrow(subset_df)

        # If fewer than 5 sites, clustering cannot proceed, groups assigned as missing
        if num_sites < 5
            println("Too few sites ($num_sites) for clustering at time step $t. Groups assigned as missing.")
            subset_df.Group = fill(missing, num_sites)
            grouping_dict[t] = subset_df
            continue  # Skip to the next time step
        end

        # Calculate distances between sites using geographical coordinates
        coordinates = select(subset_df, [:Latitude, :Longitude])
        distances = Distances.pairwise(Euclidean(), Matrix(coordinates), dims=1)
            
        # Set the initial number of clusters
        num_clusters = max(div(length(unique(subset_df.Patch)), 5), 2) # Ensure at least 2 clusters
            
        condition_met = false

        while !condition_met
            # Perform hierarchical clustering
            agglo_result = hclust(distances, linkage=:complete)
            assignments = cutree(agglo_result, k=num_clusters)
            subset_df.Group = assignments

            # Check and fix conditions for DNCI (taxa/species & sites variations)
            subset_df = Internal.check_condition_and_fix(subset_df)
            condition_met = Internal.check_conditions(subset_df)
            
            # If conditions are not met, reduce the number of clusters
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

Plots the clustering result at one time step of the `create_cluster` function using the geographic coordinates and cluster assignments of patches/sites.

Arguments
- `latitude::Vector{Float64}`: A vector of latitude coordinates of the patches/sites.
- `longitude::Vector{Float64}`: A vector of longitude coordinates of the patches/sites.
- `group::Union{AbstractVector, String}`: A vector or string indicating the cluster assignments for each data point.

Returns
- A plot showing the patches/sites colored by the cluster assignment from the `create_clusters` function.

Details
- The function assigns a unique color to each cluster and plots the patches/sites based on their geographic coordinates.
- The patches/sites are colored according to the cluster assignment.
- The plot includes black borders around the markers for better visibility.

Example
```jildoctest
julia> using MetaCommunityMetrics, Random, Distributions

julia> Random.seed!(1234)
TaskLocalRNG()

julia> time = vcat(fill(1, 20), fill(2, 20), fill(3, 20)) 
60-element Vector{Int64}:
 1
 1
 1
 ⋮
 3
 3

julia> latitude = rand(Uniform(34.0, 35.0), 60)
60-element Vector{Float64}:
 34.57986212013413
 34.41129411794985
 34.97213608245547
  ⋮
 34.57432348527832
 34.67764990759958

julia> longitude = rand(Uniform(-118.0, -117.0), 60)
60-element Vector{Float64}:
 -117.91655991056788
 -117.47420433613088
 -117.15935908052177
    ⋮
 -117.92682908601348
 -117.49824402144651

julia> patch = rand(1:10, 60) 
60-element Vector{Int64}:
 10
  2
  4
  ⋮
  6
  9

julia> total_richness = rand(5:50, 60)
60-element Vector{Int64}:
 37
 40
  6
  ⋮
 48
 28

julia> result = create_clusters(time, latitude, longitude, patch, total_richness)
Dict{Int64, DataFrame} with 3 entries:
  2 => 20×6 DataFrame…
  3 => 20×6 DataFrame…
  1 => 20×6 DataFrame…

julia> println(result[1])
20×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1   34.5799   -117.917     10              37      1
   2 │     1   34.4113   -117.474      2              40      2
   3 │     1   34.9721   -117.159      4               6      2
   4 │     1   34.0149   -117.476      9              45      1
   5 │     1   34.5204   -117.987      7              29      1
   6 │     1   34.6396   -117.594      5              18      2
   7 │     1   34.8396   -117.876     10               7      1
   8 │     1   34.9671   -117.351      1              34      2
   9 │     1   34.7898   -117.085      4              23      2
  10 │     1   34.696    -117.666      9              24      1
  11 │     1   34.5667   -117.652      8              19      1
  12 │     1   34.5364   -117.471     10              45      2
  13 │     1   34.7114   -117.042      5              43      2
  14 │     1   34.1039   -117.979      2              46      1
  15 │     1   34.8067   -117.637     10               8      2
  16 │     1   34.8705   -117.997      7              33      1
  17 │     1   34.9627   -117.369      9              48      1
  18 │     1   34.1512   -117.165      4              27      2
  19 │     1   34.7154   -117.12       8              34      1
  20 │     1   34.9395   -117.247      5              13      2

julia> result[1].Latitude
20-element Vector{Float64}:
 34.57986212013413
 34.41129411794985
 34.97213608245547
  ⋮
 34.71535467133175
 34.93954764192025

julia> result[1].Longitude
20-element Vector{Float64}:
 -117.91655991056788
 -117.47420433613088
 -117.15935908052177
    ⋮
 -117.12032610104535
 -117.24728232430559

julia> result[1].Group
20-element Vector{Int64}:
 1
 2
 2
 ⋮
 1
 2

julia> plot_clusters(result[1].Latitude, result[1].Longitude, result[1].Group)

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
    DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; count::Bool=true) -> DataFrame

Calculates the dispersal-niche continuum index (DNCI) for multiple groups, a metric proposed by Vilmi(2021). The DNCI quantifies the balance between dispersal and niche processes within a metacommunity, providing insight into community structure and the relative influence of these two key ecological drivers.

Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites or samples) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix. You can the `create_clusters` function to generate the group membership.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `count::Bool=true`: A flag indicating whether the numeber of permutations is printed. Default is `false`.

Returns
- `DataFrame`: A DataFrame containing the DNCI value, the associate confiden interval (`CI_DNCI`) and variance (`S_DNCI`) for each pair of groups.

Details
- The function calculates the DNCI for each pair of groups in the input data.
- When the DNCI value is significantly below zero, dispersal processes are likely the dominant drivers of community composition. 
- In contrast, a DNCI value significantly above zero suggests that niche processes play a primary role in shaping community composition. 
- If the DNCI value is not significantly different from zero, it indicates that dispersal and niche processes contribute equally to variations in community composition.

Example
```jildoctest
julia> using MetaCommunityMetrics, Random, Distributions

julia> Random.seed!(1234)
TaskLocalRNG()

julia> time = vcat(fill(1, 20), fill(2, 20), fill(3, 20)) 
60-element Vector{Int64}:
 1
 1
 1
 ⋮
 3
 3

julia> latitude = rand(Uniform(34.0, 35.0), 60)
60-element Vector{Float64}:
 34.57986212013413
 34.41129411794985
 34.97213608245547
  ⋮
 34.57432348527832
 34.67764990759958

julia> longitude = rand(Uniform(-118.0, -117.0), 60)
60-element Vector{Float64}:
 -117.91655991056788
 -117.47420433613088
 -117.15935908052177
    ⋮
 -117.92682908601348
 -117.49824402144651

julia> patch = rand(1:10, 60) 
60-element Vector{Int64}:
 10
  2
  4
  ⋮
  6
  9

julia> total_richness = rand(5:50, 60)
60-element Vector{Int64}:
 37
 40
  6
  ⋮
 48
 28

julia> result = create_clusters(time, latitude, longitude, patch, total_richness)
Dict{Int64, DataFrame} with 3 entries:
  2 => 20×6 DataFrame…
  3 => 20×6 DataFrame…
  1 => 20×6 DataFrame…

julia> println(result[1])
20×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1   34.5799   -117.917     10              37      1
   2 │     1   34.4113   -117.474      2              40      2
   3 │     1   34.9721   -117.159      4               6      2
   4 │     1   34.0149   -117.476      9              45      1
   5 │     1   34.5204   -117.987      7              29      1
   6 │     1   34.6396   -117.594      5              18      2
   7 │     1   34.8396   -117.876     10               7      1
   8 │     1   34.9671   -117.351      1              34      2
   9 │     1   34.7898   -117.085      4              23      2
  10 │     1   34.696    -117.666      9              24      1
  11 │     1   34.5667   -117.652      8              19      1
  12 │     1   34.5364   -117.471     10              45      2
  13 │     1   34.7114   -117.042      5              43      2
  14 │     1   34.1039   -117.979      2              46      1
  15 │     1   34.8067   -117.637     10               8      2
  16 │     1   34.8705   -117.997      7              33      1
  17 │     1   34.9627   -117.369      9              48      1
  18 │     1   34.1512   -117.165      4              27      2
  19 │     1   34.7154   -117.12       8              34      1
  20 │     1   34.9395   -117.247      5              13      2


comm = [1 0 0 1 0;
        1 1 0 0 0;
        0 1 1 0 0;
        0 0 1 1 1;
        1 0 0 0 1;
        0 1 1 0 1]

groups = ["A", "A", "B", "B", "C", "C"]
Nperm = 1000
count = false

result = DNCI_multigroup(comm, groups, Nperm; count)
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
