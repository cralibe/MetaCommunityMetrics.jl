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
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48726 rows omitted
                                                                                          
julia> total_presence_df=@pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 1, _)
791×3 DataFrame
 Row │ Species  Sampling_date_order  Total_Presence 
     │ String3  Int64                Int64          
─────┼──────────────────────────────────────────────
   1 │ BA                        41               2
   2 │ BA                        50               2
   3 │ BA                        51               8
   4 │ BA                        52              19
   5 │ BA                        53              18
  ⋮  │    ⋮              ⋮                 ⋮
 787 │ SH                        56               3
 788 │ SH                        60               4
 789 │ SH                        70               3
 790 │ SH                        73               5
 791 │ SH                       117               4
                                    781 rows omitted

julia> total_richness_df= @pipe df|>
                  innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                  groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
                  combine(_,:Presence=>sum=>:Total_Richness)|>
                  filter(row -> row[:Total_Richness] > 0, _) 

2565×5 DataFrame
  Row │ plot   Sampling_date_order  Longitude  Latitude  Total_Richness 
      │ Int64  Int64                Float64    Float64   Int64          
──────┼─────────────────────────────────────────────────────────────────
    1 │     1                   41     -110.0      35.0               5
    2 │     2                   41     -109.5      35.0               4
    3 │     4                   41     -108.5      35.0               2
    4 │     8                   41     -109.5      35.5               2
    5 │     9                   41     -109.0      35.5               3
  ⋮   │   ⋮             ⋮               ⋮         ⋮            ⋮
 2561 │     9                  117     -109.0      35.5               5
 2562 │    10                  117     -108.5      35.5               3
 2563 │    12                  117     -107.5      35.5               6
 2564 │    16                  117     -108.5      36.0               4
 2565 │    23                  117     -108.0      36.5               5
                                                       2555 rows omitted

julia> result = create_clusters(total_richness_df.Sampling_date_order, total_richness_df.Latitude, total_richness_df.Longitude, total_richness_df.plot, total_richness_df.Total_Richness)
Dict{Int64, DataFrame} with 117 entries:
  5   => 15×6 DataFrame…
  56  => 23×6 DataFrame…
  55  => 24×6 DataFrame…
  35  => 23×6 DataFrame…
  110 => 24×6 DataFrame…
  114 => 22×6 DataFrame…
  60  => 24×6 DataFrame…
  30  => 20×6 DataFrame…
  32  => 22×6 DataFrame…
  6   => 18×6 DataFrame…
  67  => 23×6 DataFrame…
  45  => 23×6 DataFrame…
  117 => 24×6 DataFrame…
  ⋮   => ⋮

julia> println(result[1])
14×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1      35.0     -110.0      1               1      1
   2 │     1      35.0     -109.5      2               1      1
   3 │     1      35.5     -109.5      8               1      1
   4 │     1      35.5     -109.0      9               1      1
   5 │     1      35.5     -108.0     11               2      2
   6 │     1      36.0     -109.5     14               2      1
   7 │     1      36.0     -108.0     17               1      2
   8 │     1      36.5     -108.5     22               1      2
   9 │     1      35.0     -107.5      6               1      2
  10 │     1      36.0     -110.0     13               1      1
  11 │     1      36.0     -109.0     15               1      1
  12 │     1      36.5     -109.5     20               1      1
  13 │     1      36.5     -109.0     21               1      2
  14 │     1      36.5     -108.0     23               1      2

```
"""
function create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}, total_richness::Vector{Int})
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Patch = patch, Total_Richness = total_richness)
    
    for t in unique(df.Time)
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
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48726 rows omitted
                                                                                          
julia> total_presence_df=@pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 1, _)
791×3 DataFrame
 Row │ Species  Sampling_date_order  Total_Presence 
     │ String3  Int64                Int64          
─────┼──────────────────────────────────────────────
   1 │ BA                        41               2
   2 │ BA                        50               2
   3 │ BA                        51               8
   4 │ BA                        52              19
   5 │ BA                        53              18
  ⋮  │    ⋮              ⋮                 ⋮
 787 │ SH                        56               3
 788 │ SH                        60               4
 789 │ SH                        70               3
 790 │ SH                        73               5
 791 │ SH                       117               4
                                    781 rows omitted

julia> total_richness_df= @pipe df|>
                  innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                  groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
                  combine(_,:Presence=>sum=>:Total_Richness)|>
                  filter(row -> row[:Total_Richness] > 0, _) 

2545×5 DataFrame
  Row │ plot   Sampling_date_order  Longitude  Latitude  Total_Richness 
      │ Int64  Int64                Float64    Float64   Int64          
──────┼─────────────────────────────────────────────────────────────────
    1 │     1                   41     -110.0      35.0               5
    2 │     2                   41     -109.5      35.0               4
    3 │     4                   41     -108.5      35.0               2
    4 │     8                   41     -109.5      35.5               2
    5 │     9                   41     -109.0      35.5               3
  ⋮   │   ⋮             ⋮               ⋮         ⋮            ⋮
 2542 │    10                  117     -108.5      35.5               3
 2543 │    12                  117     -107.5      35.5               6
 2544 │    16                  117     -108.5      36.0               4
 2545 │    23                  117     -108.0      36.5               5
                                                       2536 rows omitted

julia> clustering_result = create_clusters(total_richness_df.Sampling_date_order, total_richness_df.Latitude, total_richness_df.Longitude, total_richness_df.plot, total_richness_df.Total_Richness)
Dict{Int64, DataFrame} with 117 entries:
  5   => 17×6 DataFrame…
  56  => 23×6 DataFrame…
  55  => 24×6 DataFrame…
  35  => 23×6 DataFrame…
  110 => 24×6 DataFrame…
  114 => 22×6 DataFrame…
  60  => 24×6 DataFrame…
  30  => 20×6 DataFrame…
  32  => 22×6 DataFrame…
  6   => 19×6 DataFrame…
  67  => 23×6 DataFrame…
  45  => 23×6 DataFrame…
  117 => 24×6 DataFrame…
  73  => 23×6 DataFrame…
  ⋮   => ⋮

julia> println(clustering_result[1])
14×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1      35.0     -110.0      1               1      1
   2 │     1      35.0     -109.5      2               1      1
   3 │     1      35.5     -109.5      8               1      1
   4 │     1      35.5     -109.0      9               1      1
   5 │     1      35.5     -108.0     11               2      2
   6 │     1      36.0     -109.5     14               2      1
   7 │     1      36.0     -108.0     17               1      2
   8 │     1      36.5     -108.5     22               1      2
   9 │     1      35.0     -107.5      6               1      2
  10 │     1      36.0     -110.0     13               1      1
  11 │     1      36.0     -109.0     15               1      1
  12 │     1      36.5     -109.5     20               1      1
  13 │     1      36.5     -109.0     21               1      2
  14 │     1      36.5     -108.0     23               1      2

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

Calculates the dispersal-niche continuum index (DNCI) for multiple groups, a metric proposed by Vilmi(2021). The DNCI quantifies the balance between dispersal and niche processes within a metacommunity, providing insight into community structure and the relative influence of these two key ecological drivers. Please remove singletons (taxa/species that occuring at one patch/site within a time step) before using this function.

Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites or samples) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix. You can use the `create_clusters` function to generate the group membership.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `count::Bool=true`: A flag indicating whether the numeber of permutations is printed. Default is `false`.

Returns
- `DataFrame`: A DataFrame containing the DNCI value, the associate confiden interval (`CI_DNCI`) and variance (`S_DNCI`) for each pair of groups.

Details
- The function calculates the DNCI for each pair of groups in the input data.
- When the DNCI value is significantly below zero, dispersal processes are likely the dominant drivers of community composition. 
- In contrast, a DNCI value significantly above zero suggests that niche processes play a primary role in shaping community composition. 
- If the DNCI value is not significantly different from zero, it indicates that dispersal and niche processes contribute equally to variations in community composition.
- Please remove singletons (taxa/species that occuring at one patch/site within a time step) before using this function.
- This function is a translation/adaptation of a function from the R package `DNCImper`, licensed under GPL-3.
- Original package and documentation available at: https://github.com/Corentin-Gibert-Paleontology/DNCImper

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Random

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48726 rows omitted
                                                                                          
julia> total_presence_df=@pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 1, _)
791×3 DataFrame
 Row │ Species  Sampling_date_order  Total_Presence 
     │ String3  Int64                Int64          
─────┼──────────────────────────────────────────────
   1 │ BA                        41               2
   2 │ BA                        50               2
   3 │ BA                        51               8
   4 │ BA                        52              19
   5 │ BA                        53              18
  ⋮  │    ⋮              ⋮                 ⋮
 787 │ SH                        56               3
 788 │ SH                        60               4
 789 │ SH                        70               3
 790 │ SH                        73               5
 791 │ SH                       117               4
                                    781 rows omitted

julia> total_richness_df= @pipe df|>
                  innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                  groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
                  combine(_,:Presence=>sum=>:Total_Richness)|>
                  filter(row -> row[:Total_Richness] > 0, _) 

2565×5 DataFrame
  Row │ plot   Sampling_date_order  Longitude  Latitude  Total_Richness 
      │ Int64  Int64                Float64    Float64   Int64          
──────┼─────────────────────────────────────────────────────────────────
    1 │     1                   41     -110.0      35.0               5
    2 │     2                   41     -109.5      35.0               4
    3 │     4                   41     -108.5      35.0               2
    4 │     8                   41     -109.5      35.5               2
    5 │     9                   41     -109.0      35.5               3
  ⋮   │   ⋮             ⋮               ⋮         ⋮            ⋮
 2561 │     9                  117     -109.0      35.5               5
 2562 │    10                  117     -108.5      35.5               3
 2563 │    12                  117     -107.5      35.5               6
 2564 │    16                  117     -108.5      36.0               4
 2565 │    23                  117     -108.0      36.5               5
                                                       2555 rows omitted

julia> clustering_result = create_clusters(total_richness_df.Sampling_date_order, total_richness_df.Latitude, total_richness_df.Longitude, total_richness_df.plot, total_richness_df.Total_Richness)
Dict{Int64, DataFrame} with 117 entries:
  5   => 17×6 DataFrame…
  56  => 23×6 DataFrame…
  55  => 24×6 DataFrame…
  35  => 23×6 DataFrame…
  110 => 24×6 DataFrame…
  114 => 22×6 DataFrame…
  60  => 24×6 DataFrame…
  30  => 20×6 DataFrame…
  32  => 22×6 DataFrame…
  6   => 19×6 DataFrame…
  67  => 23×6 DataFrame…
  45  => 23×6 DataFrame…
  117 => 24×6 DataFrame…
  73  => 23×6 DataFrame…
  ⋮   => ⋮

julia> println(clustering_result[1])
14×6 DataFrame
 Row │ Time   Latitude  Longitude  Patch  Total_Richness  Group 
     │ Int64  Float64   Float64    Int64  Int64           Int64 
─────┼──────────────────────────────────────────────────────────
   1 │     1      35.0     -110.0      1               1      1
   2 │     1      35.0     -109.5      2               1      1
   3 │     1      35.5     -109.5      8               1      1
   4 │     1      35.5     -109.0      9               1      1
   5 │     1      35.5     -108.0     11               2      2
   6 │     1      36.0     -109.5     14               2      1
   7 │     1      36.0     -108.0     17               1      2
   8 │     1      36.5     -108.5     22               1      2
   9 │     1      35.0     -107.5      6               1      2
  10 │     1      36.0     -110.0     13               1      1
  11 │     1      36.0     -109.0     15               1      1
  12 │     1      36.5     -109.5     20               1      1
  13 │     1      36.5     -109.0     21               1      2
  14 │     1      36.5     -108.0     23               1      2

julia> comm= @pipe df|>
                  innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                  innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
                  filter(row -> row[:Sampling_date_order] == 1, _) |>
                  select(_, [:plot, :Species, :Presence]) |>
                  unstack(_, :Species, :Presence, fill=0) |>
                  select(_, Not(:plot)) |>
                  Matrix(_)
14×3 Matrix{Int64}:
 1  0  0
 1  0  0
 1  0  0
 1  0  0
 1  1  0
 1  1  0
 1  0  0
 0  0  1
 0  1  0
 0  1  0
 1  0  0
 0  1  0
 0  0  1
 0  0  1

julia> Random.seed!(1234) 

julia> DNCI_result = DNCI_multigroup(comm, clustering_result[1].Group, 1000; count = false)
1×5 DataFrame
 Row │ group1  group2  DNCI      CI_DNCI  S_DNCI  
     │ Int64   Int64   Float64   Float64  Float64 
─────┼────────────────────────────────────────────
   1 │      1       2  0.045603  6.52576  3.26288
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
