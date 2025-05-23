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
                    println("Warning: Cluster count fell below 2, which is not permissible for clustering. Groups assigned as missing.")
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
    plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::Union{AbstractVector, String}, output_file="clusters.svg") -> String

Visualizes clustering results by generating an SVG image displaying the geographic coordinates and cluster assignments of patches/sites.

# Arguments
- `latitude::Vector{Float64}`: A vector of latitude coordinates of the patches/sites.
- `longitude::Vector{Float64}`: A vector of longitude coordinates of the patches/sites.
- `group::Union{AbstractVector, String}`: A vector or string indicating the cluster assignments for each data point.
- `output_file::String="clusters.svg"`: The filename for the output SVG visualization. Default is "clusters.svg".

# Returns
- `String`: The path to the created SVG file.

# Details
- The function generates a standalone SVG file that can be viewed in any web browser or image viewer.
- Each cluster is assigned a unique color, and patches/sites are plotted based on their geographic coordinates.
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

julia> plot_clusters(clustering_result[1].Latitude, clustering_result[1].Longitude, clustering_result[1].Group; output_file="clusters.svg")

```
"""
function plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::Union{AbstractVector, String}; output_file="clusters.svg")
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

        GC.gc()
    end
    return ddelta
end
