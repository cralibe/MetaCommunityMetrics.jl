# src/DNCI.jl

using ..Internal


"""
    create_clusters(time::AbstractVector, latitude::Vector{Float64}, longitude::Vector{Float64}, site::AbstractVector, species::AbstractVector, presence::AbstractVector) -> Dict{Int, DataFrame}

This function creates clusters (groupings of sites) for each unique time step in a dataset which can then used for calculating DNCI. Only presnece-absence data can be used.

Arguments
- `time::AbstractVector`: Vector or single value representing sampling dates. Can be strings, integers, or any other type.
- `latitude::Vector`: A vector indicating the latitude of each site.
- `longitude::Vector`: A vector indicating the longitude of each site.
- `site::AbstractVector`: A vector indicating the spatial location of each site. At least 10 sites are required for clustering.
- `species::AbstractVector`: A vector indicating the species present at each site.
- `presence::AbstractVector`: A vector indicating the presence (1) or absence (0) of species at each site.

Returns
- `Dict{Int, DataFrame}`: A dictionary where each key represents a unique time point from the input data, with the corresponding value being a `DataFrame` for that time step. Each `DataFrame` contains the following columns: `Time`, `Latitude`, `Longitude`, `Site`, `Total_Richness`, and `Group` (indicating the assigned cluster).

Details
- This function performs hierarchical clustering on the geographical coordinates of sampling sites at each unique time step, assuming that organism dispersal occurs within the study region. 
- This function incorporates checks and adjustments to ensure the following conditions are met: at least 2 clusters, a minimum of 5 sites per cluster, and that the variation in the number of taxa/species and sites per group does not exceed 40% and 30%, respectively. These conditions are critical for calculating an unbiased DNCI value, and the function will issue warnings and the groupings will be returned as "missing" if any are not fulfilled.
- Empty sites are allowed.
- Species that is absence or presnece at all sites need to be removed before clustering, as they do not contribute to composition differences across sites.

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0                0.829467              -1.4024
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5               -1.12294               -0.0519895
     3 │  2010      1     16                    1      4  BA               0         0      35.0     -108.5               -0.409808              -0.803663
     4 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               -1.35913               -0.646369
     5 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0                0.0822                 1.09485
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 53348 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565              -0.836345
 53349 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729               -0.398522
 53350 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169               1.03257
 53351 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015               0.95971
 53352 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949               -1.59416
                                                                                                                                            53342 rows omitted
                                                                                          
julia> total_presence_df = @pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 0, _) |>
                        select(_, [:Species, :Sampling_date_order])
1038×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                        33
    2 │ BA                        41
    3 │ BA                        43
    4 │ BA                        44
    5 │ BA                        46
  ⋮   │    ⋮              ⋮
 1034 │ SH                        83
 1035 │ SH                        93
 1036 │ SH                       108
 1037 │ SH                       112
 1038 │ SH                       117
                    1028 rows omitted

julia> non_empty_site_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true)|>
                    groupby(_, [:plot]) |>
                    combine(_, :Presence=>sum=>:Total_N) |>
                    filter(row -> row[:Total_N] > 0, _) |>
                    select(_, [:plot])
24×1 DataFrame
 Row │ plot  
     │ Int64 
─────┼───────
   1 │     1
   2 │     2
   3 │     3
   4 │     4
   5 │     5
  ⋮  │   ⋮
  20 │    20
  21 │    21
  22 │    22
  23 │    23
  24 │    24
14 rows omitted

julia> ubiquitous_species_df = @pipe df |>
                        groupby(_, [:Species, :Sampling_date_order]) |>
                        combine(_, :Presence => sum => :Total_Presence) |>
                        filter(row -> row[:Total_Presence] < length(unique(non_empty_site_df.plot)), _) |>
                        select(_, [:Species, :Sampling_date_order])
2204×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                         1
    2 │ BA                         2
    3 │ BA                         3
    4 │ BA                         4
    5 │ BA                         5
  ⋮   │    ⋮              ⋮
 2200 │ SH                       113
 2201 │ SH                       114
 2202 │ SH                       115
 2203 │ SH                       116
 2204 │ SH                       117
                    2194 rows omitted

julia> filtered_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                    innerjoin(_,  non_empty_site_df, on = [:plot], makeunique = true) |>
                    innerjoin(_, ubiquitous_species_df, on = [:Species, :Sampling_date_order], makeunique = true) 

24456×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2013      4     13                   33      1  BA               0         0      35.0     -110.0                0.844714                -2.02803
     2 │  2013      4     13                   33      2  BA               0         0      35.0     -109.5                0.477563                 0.272645
     3 │  2013      4     13                   33      4  BA               0         0      35.0     -108.5                1.04654                 -0.907797
     4 │  2013      4     13                   33      8  BA               0         0      35.5     -109.5                0.170556                -0.909651
     5 │  2013      4     13                   33      9  BA               0         0      35.5     -109.0               -0.316261                -0.502008
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 24452 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565                -0.836345
 24453 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729                 -0.398522
 24454 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169                 1.03257
 24455 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015                 0.95971
 24456 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949                 -1.59416
                                                                                                                                            24446 rows omitted
julia> clustering_result = create_clusters(filtered_df.Sampling_date_order, filtered_df.Latitude, filtered_df.Longitude, filtered_df.plot, filtered_df.Species, filtered_df.Presence)
Warning: Cluster count fell below 2 at time 10, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 14, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 76, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 89, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 99, which is not permissible for clustering. Groups assigned as missing.
Dict{Int64, DataFrame} with 117 entries:
  5   => 144×7 DataFrame…
  56  => 336×7 DataFrame…
  55  => 360×7 DataFrame…
  35  => 144×7 DataFrame…
  110 => 144×7 DataFrame…
  114 => 192×7 DataFrame…
  60  => 360×7 DataFrame…
  30  => 192×7 DataFrame…
  32  => 264×7 DataFrame…
  6   => 168×7 DataFrame…
  67  => 216×7 DataFrame…
  45  => 240×7 DataFrame…
  117 => 288×7 DataFrame…
  73  => 336×7 DataFrame…
  ⋮   => ⋮

julia> clustering_result[3]
144×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │     3      35.0     -110.0      1  DM              1       1
   2 │     3      35.0     -109.5      2  DM              1       1
   3 │     3      35.0     -108.5      4  DM              0       1
   4 │     3      35.5     -109.5      8  DM              1       3
   5 │     3      35.5     -109.0      9  DM              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 140 │     3      35.5     -110.0      7  PP              1       1
 141 │     3      35.5     -108.5     10  PP              0       1
 142 │     3      36.0     -108.5     16  PP              0       2
 143 │     3      36.5     -108.0     23  PP              0       2
 144 │     3      36.5     -107.5     24  PP              0       2
                                                    134 rows omitted

julia> clustering_result[10]      
96×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group   
     │ Int64  Float64   Float64    Int64  String3  Int64     Missing 
─────┼───────────────────────────────────────────────────────────────
   1 │    10      35.0     -110.0      1  DM              1  missing 
   2 │    10      35.0     -109.5      2  DM              1  missing 
   3 │    10      35.0     -108.5      4  DM              1  missing 
   4 │    10      35.5     -109.5      8  DM              0  missing 
   5 │    10      35.5     -109.0      9  DM              1  missing 
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮         ⋮
  92 │    10      35.5     -110.0      7  OT              0  missing 
  93 │    10      35.5     -108.5     10  OT              0  missing 
  94 │    10      36.0     -108.5     16  OT              1  missing 
  95 │    10      36.5     -108.0     23  OT              0  missing 
  96 │    10      36.5     -107.5     24  OT              0  missing 
                                                      86 rows omitted                                              
```
"""
function create_clusters(time::AbstractVector, latitude::Vector{Float64}, longitude::Vector{Float64}, site::AbstractVector, species::AbstractVector, presence::AbstractVector)
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Site = site, Species = species, Presence = presence)
    
    for t in unique(df.Time)
        subset_df = filter(row -> row[:Time] == t, df)
        num_sites = length(unique(subset_df.Site))


        # If fewer than 5 sites, clustering cannot proceed, groups assigned as missing
        if num_sites < 10
            println("Too few sites ($num_sites) for clustering at time step $t. Groups assigned as missing.")
            subset_df.Group = fill(missing, nrow(subset_df))
            grouping_dict[t] = subset_df
            continue  # Skip to the next time step
        end

        # Calculate distances between sites using geographical coordinates
        coordinates = unique(select(subset_df, [:Latitude, :Longitude]))
        distances = Distances.pairwise(Euclidean(), Matrix(coordinates), dims=1)
            
        # Set the initial number of clusters
        num_clusters = max(div(length(unique(subset_df.Site)), 5), 2) # Ensure at least 2 clusters
            
        condition_met = false

        while !condition_met
            # Perform hierarchical clustering
            agglo_result = hclust(distances, linkage=:complete)
            assignments = cutree(agglo_result, k=num_clusters)
            coordinates[!, :Group] = Vector{Union{Int, Missing}}(assignments)

            # Remove any existing Group columns first
            subset_df = select(subset_df, Not(contains.(names(subset_df), "Group")))

            subset_df = @pipe subset_df |>
                innerjoin(_, coordinates, on = [:Latitude, :Longitude], makeunique = true)
        
            # Check and fix conditions with error handling
            try
                subset_df = Internal.check_condition_and_fix(subset_df)
                condition_met = Internal.check_conditions(subset_df)
            catch e
                # Log the error
                println("Error in check_condition_and_fix: ", e)
                # Set groups to missing
                subset_df.Group .= missing
                # Exit the loop by setting condition_met to true
                condition_met = true
            end
            
            # If conditions are not met, reduce the number of clusters
            if !condition_met
                num_clusters -= 1
                if num_clusters < 2
                    println("Warning: Cluster count fell below 2 at time $t, which is not permissible for clustering. Groups assigned as missing.")
                    subset_df.Group .= missing
                    # Exit the loop
                    condition_met = true
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
    plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::AbstractVector, output_file="clusters.svg") -> String

Visualizes clustering results by generating an SVG image displaying the geographic coordinates and cluster assignments of sampling sites.

# Arguments
- `latitude::Vector{Float64}`: A vector of latitude coordinates of the sampling sites.
- `longitude::Vector{Float64}`: A vector of longitude coordinates of the sampling sites.
- `group::AbstractVector`: A vector indicating the group assignments for each data point.
- `output_file::String="clusters.svg"`: The filename for the output SVG visualization. Default is "clusters.svg".

# Returns
- `String`: The path to the created SVG file.

# Details
- The function generates a standalone SVG file that can be viewed in any web browser or image viewer.
- Each cluster is assigned a unique color, and sampling sites are plotted based on their geographic coordinates.
- The visualization includes a legend identifying each cluster.


Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0                0.829467              -1.4024
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5               -1.12294               -0.0519895
     3 │  2010      1     16                    1      4  BA               0         0      35.0     -108.5               -0.409808              -0.803663
     4 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               -1.35913               -0.646369
     5 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0                0.0822                 1.09485
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 53348 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565              -0.836345
 53349 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729               -0.398522
 53350 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169               1.03257
 53351 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015               0.95971
 53352 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949               -1.59416
                                                                                                                                            53342 rows omitted
                                                                                          
julia> total_presence_df = @pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 0, _) |>
                        select(_, [:Species, :Sampling_date_order])
1038×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                        33
    2 │ BA                        41
    3 │ BA                        43
    4 │ BA                        44
    5 │ BA                        46
  ⋮   │    ⋮              ⋮
 1034 │ SH                        83
 1035 │ SH                        93
 1036 │ SH                       108
 1037 │ SH                       112
 1038 │ SH                       117
                    1028 rows omitted

julia> non_empty_site_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true)|>
                    groupby(_, [:plot]) |>
                    combine(_, :Presence=>sum=>:Total_N) |>
                    filter(row -> row[:Total_N] > 0, _) |>
                    select(_, [:plot])
24×1 DataFrame
 Row │ plot  
     │ Int64 
─────┼───────
   1 │     1
   2 │     2
   3 │     3
   4 │     4
   5 │     5
  ⋮  │   ⋮
  20 │    20
  21 │    21
  22 │    22
  23 │    23
  24 │    24
14 rows omitted

julia> ubiquitous_species_df = @pipe df |>
                        groupby(_, [:Species, :Sampling_date_order]) |>
                        combine(_, :Presence => sum => :Total_Presence) |>
                        filter(row -> row[:Total_Presence] < length(unique(non_empty_site_df.plot)), _) |>
                        select(_, [:Species, :Sampling_date_order])
2204×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                         1
    2 │ BA                         2
    3 │ BA                         3
    4 │ BA                         4
    5 │ BA                         5
  ⋮   │    ⋮              ⋮
 2200 │ SH                       113
 2201 │ SH                       114
 2202 │ SH                       115
 2203 │ SH                       116
 2204 │ SH                       117
                    2194 rows omitted

julia> filtered_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                    innerjoin(_,  non_empty_site_df, on = [:plot], makeunique = true) |>
                    innerjoin(_, ubiquitous_species_df, on = [:Species, :Sampling_date_order], makeunique = true) 

24456×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2013      4     13                   33      1  BA               0         0      35.0     -110.0                0.844714                -2.02803
     2 │  2013      4     13                   33      2  BA               0         0      35.0     -109.5                0.477563                 0.272645
     3 │  2013      4     13                   33      4  BA               0         0      35.0     -108.5                1.04654                 -0.907797
     4 │  2013      4     13                   33      8  BA               0         0      35.5     -109.5                0.170556                -0.909651
     5 │  2013      4     13                   33      9  BA               0         0      35.5     -109.0               -0.316261                -0.502008
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 24452 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565                -0.836345
 24453 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729                 -0.398522
 24454 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169                 1.03257
 24455 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015                 0.95971
 24456 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949                 -1.59416
                                                                                                                                            24446 rows omitted
julia> clustering_result = create_clusters(filtered_df.Sampling_date_order, filtered_df.Latitude, filtered_df.Longitude, filtered_df.plot, filtered_df.Species, filtered_df.Presence)
Warning: Cluster count fell below 2 at time 10, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 14, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 76, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 89, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 99, which is not permissible for clustering. Groups assigned as missing.
Dict{Int64, DataFrame} with 117 entries:
  5   => 144×7 DataFrame…
  56  => 336×7 DataFrame…
  55  => 360×7 DataFrame…
  35  => 144×7 DataFrame…
  110 => 144×7 DataFrame…
  114 => 192×7 DataFrame…
  60  => 360×7 DataFrame…
  30  => 192×7 DataFrame…
  32  => 264×7 DataFrame…
  6   => 168×7 DataFrame…
  67  => 216×7 DataFrame…
  45  => 240×7 DataFrame…
  117 => 288×7 DataFrame…
  73  => 336×7 DataFrame…
  ⋮   => ⋮

julia> clustering_result[3]
144×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │     3      35.0     -110.0      1  DM              1       1
   2 │     3      35.0     -109.5      2  DM              1       1
   3 │     3      35.0     -108.5      4  DM              0       1
   4 │     3      35.5     -109.5      8  DM              1       3
   5 │     3      35.5     -109.0      9  DM              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 140 │     3      35.5     -110.0      7  PP              1       1
 141 │     3      35.5     -108.5     10  PP              0       1
 142 │     3      36.0     -108.5     16  PP              0       2
 143 │     3      36.5     -108.0     23  PP              0       2
 144 │     3      36.5     -107.5     24  PP              0       2
                                                    134 rows omitted


julia> plot_clusters(clustering_result[3].Latitude, clustering_result[3].Longitude, clustering_result[3].Group; output_file="/Users/yc2864/Documents/research/MetaCommunityMetrics.jl/docs/src/assets/clusters.svg")

```
"""
function plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::AbstractVector; output_file="clusters.svg")
    # Get unique cluster IDs and assign numeric identifiers
    unique_clusters = unique(group)
    cluster_map = Dict(cluster => i for (i, cluster) in enumerate(unique_clusters))
    numeric_ids = [cluster_map[cluster] for cluster in group]
    
    # Define a set of distinctive colors
    base_colors = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
        "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5", 
        "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
    ]
    
    # Generate colors for each cluster
    colors = if length(unique_clusters) <= length(base_colors)
        base_colors[1:length(unique_clusters)]
    else
        # Generate additional colors if needed
        c = copy(base_colors)
        while length(c) < length(unique_clusters)
            push!(c, "#" * join(rand('0':'9', 6)))
        end
        c
    end
    
    # Calculate bounds with padding
    lon_padding = 0.05 * (maximum(longitude) - minimum(longitude))
    lat_padding = 0.05 * (maximum(latitude) - minimum(latitude))
    
    min_lon = minimum(longitude) - lon_padding
    max_lon = maximum(longitude) + lon_padding
    min_lat = minimum(latitude) - lat_padding
    max_lat = maximum(latitude) + lat_padding
    
    # SVG dimensions
    width = 800
    height = 600
    margin = 100
    
    # Map coordinates to SVG space
    function map_coords(lon, lat)
        x = margin + (lon - min_lon) / (max_lon - min_lon) * (width - 2 * margin)
        # Flip y-axis (SVG has origin at top-left)
        y = height - margin - (lat - min_lat) / (max_lat - min_lat) * (height - 2 * margin)
        return round(x, digits=1), round(y, digits=1)
    end
    
    # Start building SVG content
    svg = """<?xml version="1.0" encoding="UTF-8"?>
    <svg width="$(width)" height="$(height)" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    
    <!-- Title -->
    <text x="$(width/2)" y="25" font-family="Arial" font-size="18" text-anchor="middle" font-weight="bold">Cluster Visualization</text>
    
    <!-- Axes -->
    <line x1="$(margin)" y1="$(height-margin)" x2="$(width-margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    <line x1="$(margin)" y1="$(margin)" x2="$(margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    
    <!-- Axis labels -->
    <text x="$(width/2)" y="$(height-10)" font-family="Arial" font-size="14" text-anchor="middle">Longitude</text>
    <text x="15" y="$(height/2)" font-family="Arial" font-size="14" text-anchor="middle" transform="rotate(-90, 15, $(height/2))">Latitude</text>
    
    <!-- Data points -->
    """
    
    # Add all data points
    for i in 1:length(latitude)
        x, y = map_coords(longitude[i], latitude[i])
        cluster_id = numeric_ids[i]
        color = colors[cluster_id]
        
        svg *= """
        <circle cx="$(x)" cy="$(y)" r="5" fill="$(color)" stroke="black" stroke-width="0.5" opacity="0.8"/>
        """
    end
    
    # Add legend
    legend_x = width - margin + 15
    legend_y = margin + 20
    
    svg *= """
    <!-- Legend -->
    <text x="$(legend_x)" y="$(legend_y - 15)" font-family="Arial" font-size="12" font-weight="bold">Clusters</text>
    """
    
    for (i, cluster) in enumerate(unique_clusters)
        y_pos = legend_y + (i-1) * 20
        svg *= """
        <rect x="$(legend_x)" y="$(y_pos-10)" width="10" height="10" fill="$(colors[i])" stroke="black" stroke-width="0.5"/>
        <text x="$(legend_x + 20)" y="$(y_pos)" font-family="Arial" font-size="12">$(cluster)</text>
        """
    end
    
    # Close SVG
    svg *= "</svg>"
    
    # Write to file
    open(output_file, "w") do f
        write(f, svg)
    end
    
    println("Cluster visualization saved to $(abspath(output_file))")
    println("Open this file in any web browser or image viewer to see the visualization")
    
    return output_file
end


"""
    DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; Nperm_count::Bool=true) -> DataFrame

Calculates the dispersal-niche continuum index (DNCI) for multiple groups, a metric proposed by Vilmi(2021). The DNCI quantifies the balance between dispersal and niche processes within a metacommunity, providing insight into community structure and the relative influence of these two key ecological drivers. 

Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites or samples) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix. You can use the `create_clusters` function to generate the group membership.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `Nperm_count::Bool=true`: A flag indicating whether the number of permutations is printed. Default is `false`.

Returns
- The DataFrame will have the following columns:
  - `Group1`: The first group in the pair.
  - `Group2`: The second group in the pair.
  - `DNCI`: The calculated DNCI value.
  - `CI_DNCI`: The confidence interval for the DNCI value.
  - `S_DNCI`: The variance of the DNCI value.
  - `Status`: A string indicating whether how the DNCI is calculated. It is mainly used to flag edge cases.
    - `normal` indicates that the DNCI is calculated as normal.
    - `empty_community` indicates no species existed in any sites in a given group pair, `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `only_one_species_exists` indicates that only one species existed in a given group pair, which is not possible to calculate relative species contribution to overall dissimilarity. `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `quasi_swap_permutation_not_possible` indicates that the quasi-swap permutation (a matrix permutation algorithms that preserves row and column sums) is not possible due to extreme matrix constraints that prevent any rearrangement of species across sites. DNCI, CI_DNCI, and S_DNCI are returned as NaN.
    - `one_way_to_quasi_swap` indicates that only one arrangement is possible under quasi-swap constraints, preventing generation of a null distribution. DNCI, CI_DNCI, and S_DNCI are returned as NaN.
    - `inadequate_variation_quasi_swap` indicates that quasi-swap permutations generated insufficient variation (coefficient of variation <1%) for reliable statistical inference. DNCI, CI_DNCI, and S_DNCI are returned as NaN.
Details
- The function calculates the DNCI for each pair of groups in the input data.
- When the DNCI value is significantly below zero, dispersal processes are likely the dominant drivers of community composition. 
- In contrast, a DNCI value significantly above zero suggests that niche processes play a primary role in shaping community composition. 
- If the DNCI value is not significantly different from zero, it indicates that dispersal and niche processes contribute equally to spatial variations in community composition at a given time point.
- Different from the original implementation, empty sites are allowed. However, species that is absence or presnece at all sites need to be removed before clustering, as they do not contribute to composition differences across sites.
- This function is a adaptation of the function `DNCI_multigroup()` from the R package `DNCImper`, licensed under GPL-3.
- Original package and documentation available at: https://github.com/Corentin-Gibert-Paleontology/DNCImper


Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Random

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0                0.829467              -1.4024
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5               -1.12294               -0.0519895
     3 │  2010      1     16                    1      4  BA               0         0      35.0     -108.5               -0.409808              -0.803663
     4 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               -1.35913               -0.646369
     5 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0                0.0822                 1.09485
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 53348 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565              -0.836345
 53349 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729               -0.398522
 53350 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169               1.03257
 53351 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015               0.95971
 53352 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949               -1.59416
                                                                                                                                            53342 rows omitted
                                                                                          
julia> total_presence_df = @pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 0, _) |>
                        select(_, [:Species, :Sampling_date_order])
1038×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                        33
    2 │ BA                        41
    3 │ BA                        43
    4 │ BA                        44
    5 │ BA                        46
  ⋮   │    ⋮              ⋮
 1034 │ SH                        83
 1035 │ SH                        93
 1036 │ SH                       108
 1037 │ SH                       112
 1038 │ SH                       117
                    1028 rows omitted

julia> non_empty_site_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true)|>
                    groupby(_, [:plot]) |>
                    combine(_, :Presence=>sum=>:Total_N) |>
                    filter(row -> row[:Total_N] > 0, _) |>
                    select(_, [:plot])
24×1 DataFrame
 Row │ plot  
     │ Int64 
─────┼───────
   1 │     1
   2 │     2
   3 │     3
   4 │     4
   5 │     5
  ⋮  │   ⋮
  20 │    20
  21 │    21
  22 │    22
  23 │    23
  24 │    24
14 rows omitted

julia> ubiquitous_species_df = @pipe df |>
                        groupby(_, [:Species, :Sampling_date_order]) |>
                        combine(_, :Presence => sum => :Total_Presence) |>
                        filter(row -> row[:Total_Presence] < length(unique(non_empty_site_df.plot)), _) |>
                        select(_, [:Species, :Sampling_date_order])
2204×2 DataFrame
  Row │ Species  Sampling_date_order 
      │ String3  Int64               
──────┼──────────────────────────────
    1 │ BA                         1
    2 │ BA                         2
    3 │ BA                         3
    4 │ BA                         4
    5 │ BA                         5
  ⋮   │    ⋮              ⋮
 2200 │ SH                       113
 2201 │ SH                       114
 2202 │ SH                       115
 2203 │ SH                       116
 2204 │ SH                       117
                    2194 rows omitted

julia> filtered_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                    innerjoin(_,  non_empty_site_df, on = [:plot], makeunique = true) |>
                    innerjoin(_, ubiquitous_species_df, on = [:Species, :Sampling_date_order], makeunique = true) 

24456×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2013      4     13                   33      1  BA               0         0      35.0     -110.0                0.844714                -2.02803
     2 │  2013      4     13                   33      2  BA               0         0      35.0     -109.5                0.477563                 0.272645
     3 │  2013      4     13                   33      4  BA               0         0      35.0     -108.5                1.04654                 -0.907797
     4 │  2013      4     13                   33      8  BA               0         0      35.5     -109.5                0.170556                -0.909651
     5 │  2013      4     13                   33      9  BA               0         0      35.5     -109.0               -0.316261                -0.502008
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 24452 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565                -0.836345
 24453 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729                 -0.398522
 24454 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169                 1.03257
 24455 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015                 0.95971
 24456 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949                 -1.59416
                                                                                                                                            24446 rows omitted
julia> clustering_result = create_clusters(filtered_df.Sampling_date_order, filtered_df.Latitude, filtered_df.Longitude, filtered_df.plot, filtered_df.Species, filtered_df.Presence)
Warning: Cluster count fell below 2 at time 10, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 14, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 76, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 89, which is not permissible for clustering. Groups assigned as missing.
Warning: Cluster count fell below 2 at time 99, which is not permissible for clustering. Groups assigned as missing.
Dict{Int64, DataFrame} with 117 entries:
  5   => 144×7 DataFrame…
  56  => 336×7 DataFrame…
  55  => 360×7 DataFrame…
  35  => 144×7 DataFrame…
  110 => 144×7 DataFrame…
  114 => 192×7 DataFrame…
  60  => 360×7 DataFrame…
  30  => 192×7 DataFrame…
  32  => 264×7 DataFrame…
  6   => 168×7 DataFrame…
  67  => 216×7 DataFrame…
  45  => 240×7 DataFrame…
  117 => 288×7 DataFrame…
  73  => 336×7 DataFrame…
  ⋮   => ⋮

julia> clustering_result[3]
144×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │     3      35.0     -110.0      1  DM              1       1
   2 │     3      35.0     -109.5      2  DM              1       1
   3 │     3      35.0     -108.5      4  DM              0       1
   4 │     3      35.5     -109.5      8  DM              1       3
   5 │     3      35.5     -109.0      9  DM              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 140 │     3      35.5     -110.0      7  PP              1       1
 141 │     3      35.5     -108.5     10  PP              0       1
 142 │     3      36.0     -108.5     16  PP              0       2
 143 │     3      36.5     -108.0     23  PP              0       2
 144 │     3      36.5     -107.5     24  PP              0       2
                                                    134 rows omitted

julia> group_df = @pipe filtered_df |>
                  filter(row -> row[:Sampling_date_order] == 3, _) |>
                  select(_, [:plot, :Species, :Presence]) |>
                  innerjoin(_, clustering_result[3], on = [:plot => :Site, :Species], makeunique = true)|>
                  select(_, [:plot, :Species, :Presence, :Group]) |>
                  unstack(_, :Species, :Presence, fill=0)
24×8 DataFrame
 Row │ plot   Group   DM     DO     OL     OT     PB     PP    
     │ Int64  Int64?  Int64  Int64  Int64  Int64  Int64  Int64 
─────┼─────────────────────────────────────────────────────────
   1 │     1       1      1      0      0      0      0      0
   2 │     2       1      1      0      0      0      1      0
   3 │     4       1      0      0      0      0      0      1
   4 │     8       3      1      1      0      0      0      0
   5 │     9       1      0      0      0      0      0      0
  ⋮  │   ⋮      ⋮       ⋮      ⋮      ⋮      ⋮      ⋮      ⋮
  20 │     7       1      0      0      0      1      0      1
  21 │    10       1      0      0      0      0      0      0
  22 │    16       2      0      0      0      0      0      0
  23 │    23       2      0      0      0      1      0      0
  24 │    24       2      0      0      0      0      0      0
                                                14 rows omitted
                                                                                                                                           5 rows omitted

julia> comm= @pipe group_df |>
                  select(_, Not([:plot,:Group])) |>
                  Matrix(_)
24×6 Matrix{Int64}:
 1  0  0  0  0  0
 1  0  0  0  1  0
 0  0  0  0  0  1
 1  1  0  0  0  0
 0  0  0  0  0  0
 1  0  0  0  0  0
 0  0  0  0  0  0
 ⋮              ⋮
 0  0  0  0  1  0
 0  0  0  0  0  0
 0  0  0  1  0  1
 0  0  0  0  0  0
 0  0  0  0  0  0
 0  0  0  1  0  0
 0  0  0  0  0  0

julia> Random.seed!(1234) 

julia> DNCI_result = DNCI_multigroup(comm, group_df.Group, 1000; Nperm_count = false)
3×6 DataFrame
 Row │ group1  group2  DNCI      CI_DNCI  S_DNCI   status 
     │ Int64   Int64   Float64   Float64  Float64  String 
─────┼────────────────────────────────────────────────────
   1 │      1       2  -1.84521  3.22534  1.61267  normal
   2 │      1       3  -3.47866  2.33283  1.16642  normal
   3 │      2       3  -2.85858  2.28596  1.14298  normal
```
"""
function DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; Nperm_count::Bool=true) #for presence-absence data only
    
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

        # Calculate logical arrays for filtering
        non_zero_sum_columns = sum(paired_x, dims=1) .!= 0
        non_ubiquitous_columns = sum(paired_x, dims=1) .!= size(paired_x, 1)  # Not present everywhere

        # Combine both conditions: non-zero AND not ubiquitous
        valid_columns = non_zero_sum_columns .& non_ubiquitous_columns

        # Convert logical index to actual column indices
        column_indices = findall(x -> x, valid_columns[:])

        # Check if any columns remain
        if isempty(column_indices)
            println("Empty community detected, returning a metric with NaN values.")
            DNCI_result = DataFrame(group1=group_combinations[i][1], group2=group_combinations[i][2], 
                        DNCI=NaN, CI_DNCI=NaN, S_DNCI=NaN,
                        status="empty_community")
            return DNCI_result
        end
    
        # Subset the matrix using these indices
        paired_x = paired_x[:, column_indices]

        # Ensure that the paired_x matrix has at least two species
        if size(paired_x, 2) == 1
            println("Only one species present, not enough species in the selected groups to calculate DNCI. Skipping this group pair.")
            DNCI_result = DataFrame(group1=group_combinations[i][1], group2=group_combinations[i][2], 
                          DNCI=NaN, CI_DNCI=NaN, S_DNCI=NaN,
                          status= "only_one_species_exists")
            return DNCI_result
        end


        group_pair = vcat(
            fill(group_combinations[i][1], size(splitx[group_combinations[i][1]], 1)),
            fill(group_combinations[i][2], size(splitx[group_combinations[i][2]], 1)))

        DNCI_result = Internal.DNCI_ses(paired_x, group_pair, Nperm; count=Nperm_count)
            
        append!(ddelta, DNCI_result, promote=true)

        GC.gc()
    end
    return ddelta
end
