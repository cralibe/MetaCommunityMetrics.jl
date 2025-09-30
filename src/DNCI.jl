# src/DNCI.jl

using ..Internal


"""
    create_groups(time::AbstractVector, latitude::Vector{Float64}, longitude::Vector{Float64}, site::AbstractVector, species::AbstractVector, presence::AbstractVector) -> Dict{Int, DataFrame}

This function creates groupings of sites for each unique time step in a dataset which can then used for calculating DNCI. Only presnece-absence data can be used.

Arguments
- `time::AbstractVector`: Vector or single value representing sampling dates. Can be strings, integers, or any other type.
- `latitude::Vector`: A vector indicating the latitude of each site.
- `longitude::Vector`: A vector indicating the longitude of each site.
- `site::AbstractVector`: A vector indicating the spatial location of each site. At least 10 sites are required for clustering.
- `species::AbstractVector`: A vector indicating the species present at each site.
- `presence::AbstractVector`: A vector indicating the presence (1) or absence (0) of species at each site.

Returns
- `Dict{Int, DataFrame}`: A dictionary where each key represents a unique time point from the input data, with the corresponding value being a `DataFrame` for that time step. Each `DataFrame` contains the following columns:
    - `Time`
    - `Latitude`
    - `Longitude`
    - `Site`
    - `Species`
    - `Presence`
    - `Group` (indicating the assigned groups).

Details
- This function performs hierarchical clustering on the geographical coordinates of sampling sites at each time point separately and processes all time points in a single execution. 
- This function incorporates checks and adjustments to ensure the following conditions are met: 
    - Having at least 2 groups
    - A minimum of 5 sites per group, 
    - The variation in the number of taxa/species and sites per group does not exceed 40% and 30%, respectively. 
    These conditions are critical for calculating an unbiased DNCI value, and the function will issue warnings and the groupings will be returned as "missing" if any of the above are not fulfilled.
- Empty sites are allowed.
Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
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
                                                                                          
julia> grouping_result = create_groups(df.Sampling_date_order, df.Latitude, df.Longitude, df.plot, df.Species, df.Presence)
Warning: Group count fell below 2 at time 10, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 14, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 76, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 89, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 99, which is not permissible for DNCI analysis. Groups assigned as missing.
Dict{Int64, DataFrames.DataFrame} with 117 entries:
  5   => 456×7 DataFrame…
  56  => 456×7 DataFrame…
  35  => 456×7 DataFrame…
  55  => 456×7 DataFrame…
  110 => 456×7 DataFrame…
  114 => 456×7 DataFrame…
  60  => 456×7 DataFrame…
  30  => 456×7 DataFrame…
  32  => 456×7 DataFrame…
  6   => 456×7 DataFrame…
  67  => 456×7 DataFrame…
  45  => 456×7 DataFrame…
  117 => 456×7 DataFrame…
  73  => 456×7 DataFrame…
  ⋮   => ⋮

julia> grouping_result[10]      
456×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group   
     │ Int64  Float64   Float64    Int64  String3  Int64     Missing 
─────┼───────────────────────────────────────────────────────────────
   1 │    10      35.0     -110.0      1  BA              0  missing 
   2 │    10      35.0     -109.5      2  BA              0  missing 
   3 │    10      35.0     -108.5      4  BA              0  missing 
   4 │    10      35.5     -109.5      8  BA              0  missing 
   5 │    10      35.5     -109.0      9  BA              0  missing 
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮         ⋮
 452 │    10      35.5     -110.0      7  SH              0  missing 
 453 │    10      35.5     -108.5     10  SH              0  missing 
 454 │    10      36.0     -108.5     16  SH              0  missing 
 455 │    10      36.5     -108.0     23  SH              0  missing 
 456 │    10      36.5     -107.5     24  SH              0  missing 
                                                     446 rows omitted      
 
julia> grouping_result[60]                                                      
456×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │    60      35.0     -108.5      4  BA              0       1
   2 │    60      35.0     -108.0      5  BA              1       1
   3 │    60      35.0     -107.5      6  BA              0       1
   4 │    60      35.5     -110.0      7  BA              0       2
   5 │    60      35.5     -108.0     11  BA              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 452 │    60      35.5     -109.0      9  SH              0       2
 453 │    60      35.5     -108.5     10  SH              0       1
 454 │    60      35.5     -107.5     12  SH              1       1
 455 │    60      36.0     -108.5     16  SH              0       4
 456 │    60      36.5     -108.0     23  SH              0       4
                                                    446 rows omitted
```
"""
function create_groups(time::AbstractVector, latitude::Vector{Float64}, longitude::Vector{Float64}, site::AbstractVector, species::AbstractVector, presence::AbstractVector)
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Site = site, Species = species, Presence = presence)
    
    for t in unique(df.Time)
        subset_df = filter(row -> row[:Time] == t, df)
        num_sites = length(unique(subset_df.Site))


        # If fewer than 5 sites, clustering cannot proceed, groups assigned as missing
        if num_sites < 10
            println("Too few sites ($num_sites) for grouping at time step $t. Groups assigned as missing.")
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
                    println("Warning: Group count fell below 2 at time $t, which is not permissible for DNCI analysis. Groups assigned as missing.")
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
    plot_groups(latitude::Vector{Float64}, longitude::Vector{Float64}, group::AbstractVector, output_file="groups.svg") -> String

Visualizes grouping results by generating an SVG image displaying the geographic coordinates and cluster assignments of sampling sites.

# Arguments
- `latitude::Vector{Float64}`: A vector of latitude coordinates of the sampling sites.
- `longitude::Vector{Float64}`: A vector of longitude coordinates of the sampling sites.
- `group::AbstractVector`: A vector indicating the group assignments for each data point.
- `output_file::String="clusters.svg"`: The filename for the output SVG visualization. Default is "groups.svg".

# Returns
- `String`: The path to the created SVG file.

# Details
- The functions provides visualization for one time point per function call.
- The function generates a standalone SVG file that can be viewed in any web browser or image viewer.
- Each group is assigned a unique color, and sampling sites are plotted based on their geographic coordinates.
- The visualization includes a legend identifying each group.


Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
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
                                                                                          
julia> grouping_result = create_groups(df.Sampling_date_order, df.Latitude, df.Longitude, df.plot, df.Species, df.Presence)
Warning: Group count fell below 2 at time 10, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 14, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 76, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 89, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 99, which is not permissible for DNCI analysis. Groups assigned as missing.
Dict{Int64, DataFrames.DataFrame} with 117 entries:
  5   => 456×7 DataFrame…
  56  => 456×7 DataFrame…
  35  => 456×7 DataFrame…
  55  => 456×7 DataFrame…
  110 => 456×7 DataFrame…
  114 => 456×7 DataFrame…
  60  => 456×7 DataFrame…
  30  => 456×7 DataFrame…
  32  => 456×7 DataFrame…
  6   => 456×7 DataFrame…
  67  => 456×7 DataFrame…
  45  => 456×7 DataFrame…
  117 => 456×7 DataFrame…
  73  => 456×7 DataFrame…
  ⋮   => ⋮

julia> grouping_result[60]                                                      
456×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │    60      35.0     -108.5      4  BA              0       1
   2 │    60      35.0     -108.0      5  BA              1       1
   3 │    60      35.0     -107.5      6  BA              0       1
   4 │    60      35.5     -110.0      7  BA              0       2
   5 │    60      35.5     -108.0     11  BA              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 452 │    60      35.5     -109.0      9  SH              0       2
 453 │    60      35.5     -108.5     10  SH              0       1
 454 │    60      35.5     -107.5     12  SH              1       1
 455 │    60      36.0     -108.5     16  SH              0       4
 456 │    60      36.5     -108.0     23  SH              0       4
                                                    446 rows omitted

julia> plot_groups(grouping_result[60].Latitude, grouping_result[60].Longitude, grouping_result[60].Group; output_file="groups.svg")

```
"""
function plot_groups(latitude::Vector{Float64}, longitude::Vector{Float64}, group::AbstractVector; output_file="groups.svg")
    # Get unique cluster IDs and assign numeric identifiers
    unique_clusters = unique(group)
    cluster_map = Dict(cluster => i for (i, cluster) in enumerate(unique_clusters))
    numeric_ids = [cluster_map[cluster] for cluster in group]
    
    # define a set of distinctive colors
    base_colors = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
        "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5", 
        "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
    ]
    
    # generate colors for each cluster
    colors = if length(unique_clusters) <= length(base_colors)
        base_colors[1:length(unique_clusters)]
    else
        # G\generate additional colors if needed
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
    
    # image dimensions
    width = 800
    height = 600
    margin = 100
    
    # map coordinates to SVG space
    function map_coords(lon, lat)
        x = margin + (lon - min_lon) / (max_lon - min_lon) * (width - 2 * margin)
        # Flip y-axis (SVG has origin at top-left)
        y = height - margin - (lat - min_lat) / (max_lat - min_lat) * (height - 2 * margin)
        return round(x, digits=1), round(y, digits=1)
    end
    
    # the SVG content
    svg = """<?xml version="1.0" encoding="UTF-8"?>
    <svg width="$(width)" height="$(height)" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    
    <!-- Title -->
    <text x="$(width/2)" y="25" font-family="Arial" font-size="20" text-anchor="middle" font-weight="bold">Grouping Result</text>
    
    <!-- Axes -->
    <line x1="$(margin)" y1="$(height-margin)" x2="$(width-margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    <line x1="$(margin)" y1="$(margin)" x2="$(margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    
    <!-- Axis labels -->
    <text x="$(width/2)" y="$(height-10)" font-family="Arial" font-size="18" text-anchor="middle">Longitude</text>
    <text x="15" y="$(height/2)" font-family="Arial" font-size="16" text-anchor="middle" transform="rotate(-90, 15, $(height/2))">Latitude</text>
    
    <!-- Data points -->
    """
    
    # add all data points
    for i in 1:length(latitude)
        x, y = map_coords(longitude[i], latitude[i])
        cluster_id = numeric_ids[i]
        color = colors[cluster_id]
        
        svg *= """
        <circle cx="$(x)" cy="$(y)" r="5" fill="$(color)" stroke="black" stroke-width="0.5" opacity="0.8"/>
        """
    end
    
    # add legend
    legend_x = width - margin + 15
    legend_y = margin + 20
    
    svg *= """
    <!-- Legend -->
    <text x="$(legend_x)" y="$(legend_y - 15)" font-family="Arial" font-size="18" font-weight="bold">Groups</text>
    """
    
    for (i, cluster) in enumerate(unique_clusters)
        y_pos = legend_y + (i-1) * 20
        svg *= """
        <rect x="$(legend_x)" y="$(y_pos-10)" width="10" height="10" fill="$(colors[i])" stroke="black" stroke-width="0.5"/>
        <text x="$(legend_x + 20)" y="$(y_pos)" font-family="Arial" font-size="16">$(cluster)</text>
        """
    end
    
    # close SVG
    svg *= "</svg>"
    
    # write to file
    open(output_file, "w") do f
        write(f, svg)
    end
    
    println("Cluster visualization saved to $(abspath(output_file))")
    println("Open this file in any web browser or image viewer to see the visualization")
    
    return output_file
end


"""
    DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; Nperm_count::Bool=true) -> DataFrame

Calculates the dispersal-niche continuum index (DNCI) for a metacommunity, a metric proposed by Vilmi(2021). The DNCI quantifies the balance between dispersal and niche processes within a metacommunity, providing insight into community structure and the relative influence of these two key ecological drivers. 

Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix. You can use the `create_clusters` function to generate the group membership.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `Nperm_count::Bool=true`: A flag indicating whether the number of permutations is printed. Default is `false`.

Returns
The DataFrame will have the following columns:
- `Group1`: The first group in the pair.
- `Group2`: The second group in the pair.
- `DNCI`: The calculated DNCI value.
- `CI_DNCI`: The confidence interval for the DNCI value.
- `S_DNCI`: The standard deviation of the DNCI value.
- `Status`: A string indicating how the DNCI is calculated. It is mainly used to flag edge cases as follows:
    - `normal` indicates that the DNCI is calculated as normal.
    - `empty_community` indicates no species existed in any sites in a given group pair, `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `only_one_species_exists` indicates that only one species existed in a given group pair, which is not possible to calculate relative species contribution to overall dissimilarity. `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `quasi_swap_permutation_not_possible` indicates that the quasi-swap permutation (a matrix permutation algorithms that preserves row and column sums) is not possible due to extreme matrix constraints that prevent any rearrangement of species across sites. `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `one_way_to_quasi_swap` indicates that only one arrangement is possible under quasi-swap constraints, preventing generation of a null distribution. `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
    - `inadequate_variation_quasi_swap` indicates that quasi-swap permutations generated insufficient variation (coefficient of variation <1%) for reliable statistical inference. `DNCI`, `CI_DNCI`, and `S_DNCI` are returned as `NaN`.
Details
- The function calculates the DNCI for each pair of groups in the input data.
- When the DNCI value is significantly below zero, dispersal processes are likely the dominant drivers of community composition. 
- In contrast, a DNCI value significantly above zero suggests that niche processes play a primary role in shaping community composition. 
- If the DNCI value is not significantly different from zero, it indicates that dispersal and niche processes contribute equally to spatial variations in community composition at a given time point.
- Different from the original implementation, empty sites and singletons (species that only occupy one site at a given time) are allowed.
- This function is a adaptation of the function `DNCI_multigroup()` from the R package `DNCImper`, licensed under GPL-3.
- Original package and documentation available at: https://github.com/Corentin-Gibert-Paleontology/DNCImper


Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Random

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
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
                                                                                          
julia> grouping_result = create_groups(df.Sampling_date_order, df.Latitude, df.Longitude, df.plot, df.Species, df.Presence)
Warning: Group count fell below 2 at time 10, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 14, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 76, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 89, which is not permissible for DNCI analysis. Groups assigned as missing.
Warning: Group count fell below 2 at time 99, which is not permissible for DNCI analysis. Groups assigned as missing.
Dict{Int64, DataFrames.DataFrame} with 117 entries:
  5   => 456×7 DataFrame…
  56  => 456×7 DataFrame…
  35  => 456×7 DataFrame…
  55  => 456×7 DataFrame…
  110 => 456×7 DataFrame…
  114 => 456×7 DataFrame…
  60  => 456×7 DataFrame…
  30  => 456×7 DataFrame…
  32  => 456×7 DataFrame…
  6   => 456×7 DataFrame…
  67  => 456×7 DataFrame…
  45  => 456×7 DataFrame…
  117 => 456×7 DataFrame…
  73  => 456×7 DataFrame…
  ⋮   => ⋮

julia> grouping_result[60]
456×7 DataFrame
 Row │ Time   Latitude  Longitude  Site   Species  Presence  Group  
     │ Int64  Float64   Float64    Int64  String3  Int64     Int64? 
─────┼──────────────────────────────────────────────────────────────
   1 │    60      35.0     -108.5      4  BA              0       1
   2 │    60      35.0     -108.0      5  BA              1       1
   3 │    60      35.0     -107.5      6  BA              0       1
   4 │    60      35.5     -110.0      7  BA              0       2
   5 │    60      35.5     -108.0     11  BA              0       1
  ⋮  │   ⋮       ⋮          ⋮        ⋮       ⋮        ⋮        ⋮
 452 │    60      35.5     -109.0      9  SH              0       2
 453 │    60      35.5     -108.5     10  SH              0       1
 454 │    60      35.5     -107.5     12  SH              1       1
 455 │    60      36.0     -108.5     16  SH              0       4
 456 │    60      36.5     -108.0     23  SH              0       4
                                                    446 rows omitted

julia> group_df = @pipe df |>
                filter(row -> row[:Sampling_date_order] == 60, _) |>
                select(_, [:plot, :Species, :Presence]) |>
                innerjoin(_, grouping_result[60], on = [:plot => :Site, :Species], makeunique = true)|>
                select(_, [:plot, :Species, :Presence, :Group]) |>
                unstack(_, :Species, :Presence, fill=0)
24×21 DataFrame
 Row │ plot   Group   BA     DM     DO     DS     NA     OL     OT     PB     PE     PF     PH     PL     PM     PP     RF     RM     RO     SF     SH    
     │ Int64  Int64?  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64 
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │     4       1      0      1      1      0      0      1      1      0      0      0      0      0      0      0      0      1      0      1      0
   2 │     5       1      1      1      1      0      0      0      1      0      1      0      0      0      0      1      0      1      0      0      0
   3 │     6       1      0      1      1      0      0      0      0      0      1      0      0      0      0      1      0      1      0      0      0
   4 │     7       2      0      1      1      0      0      1      1      0      0      0      0      0      0      1      0      0      0      0      1
   5 │    11       1      0      1      1      0      0      0      1      0      0      0      0      0      0      1      0      0      0      0      0
  ⋮  │   ⋮      ⋮       ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮      ⋮
  20 │     9       2      0      1      0      0      0      1      1      0      1      0      0      0      0      0      0      0      0      0      0
  21 │    10       1      0      0      0      0      0      0      0      0      1      0      0      0      0      1      0      0      0      0      0
  22 │    12       1      0      0      0      0      1      0      1      0      1      0      0      0      1      0      0      1      0      0      1
  23 │    16       4      0      0      1      0      0      0      1      0      1      0      0      0      0      0      0      1      0      0      0
  24 │    23       4      0      1      0      0      0      0      0      1      0      0      0      0      1      0      0      0      0      0      0
                                                                                                                                           14 rows omitted
                                                                                                                                          
julia> comm= @pipe group_df |>
                  select(_, Not([:plot,:Group])) |>
                  Matrix(_)
24×19 Matrix{Int64}:
 0  1  1  0  0  1  1  0  0  0  0  0  0  0  0  1  0  1  0
 1  1  1  0  0  0  1  0  1  0  0  0  0  1  0  1  0  0  0
 0  1  1  0  0  0  0  0  1  0  0  0  0  1  0  1  0  0  0
 0  1  1  0  0  1  1  0  0  0  0  0  0  1  0  0  0  0  1
 0  1  1  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0
 0  1  1  0  0  0  0  1  1  0  0  0  0  0  0  1  0  0  0
 1  1  1  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0
 ⋮              ⋮              ⋮              ⋮        
 0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  1  0  0  1  1  0  0  0  0  0  0
 0  1  0  0  0  1  1  0  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0  0
 0  0  0  0  1  0  1  0  1  0  0  0  1  0  0  1  0  0  1
 0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  1  0  0  0
 0  1  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0  0  0

julia> Random.seed!(1234) 

julia> DNCI_result = DNCI_multigroup(comm, group_df.Group, 1000; Nperm_count = false)
6×6 DataFrame
 Row │ group1  group2  DNCI      CI_DNCI  S_DNCI    status 
     │ Int64   Int64   Float64   Float64  Float64   String 
─────┼─────────────────────────────────────────────────────
   1 │      1       2  -3.41127  2.17348  1.08674   normal
   2 │      1       3  -2.44866  2.05951  1.02976   normal
   3 │      1       4  -2.3671   2.45697  1.22848   normal
   4 │      2       3  -2.65022  2.28931  1.14466   normal
   5 │      2       4  -3.0168   2.43496  1.21748   normal
   6 │      3       4  -1.83521  1.9589   0.979449  normal
```
"""
function DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000; Nperm_count::Bool=true) #for presence-absence data only
    
    group_combinations = collect(combinations(unique(sort(groups)),2))

    ddelta = DataFrame()

    for i in 1:size(group_combinations,1)
        # Create an empty dictionary to hold the split data
        splitx = Dict()

        # comm is a matrix and groups is a vector indicating the group for each row in comm
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

        # access and concatenate matrices from two groups
        if haskey(splitx, group_combinations[i][1]) && haskey(splitx, group_combinations[i][2])
            paired_x = vcat(splitx[group_combinations[i][1]], splitx[group_combinations[i][2]])
        else
            error("One of the groups does not exist in the dictionary.")
        end

        # Calculate logical arrays for filtering
        non_zero_sum_columns = sum(paired_x, dims=1) .!= 0
        non_ubiquitous_columns = sum(paired_x, dims=1) .!= size(paired_x, 1)  # Not present everywhere

        # Combine both conditions: non-zero and not ubiquitous
        valid_columns = non_zero_sum_columns .& non_ubiquitous_columns

        # Convert logical index to actual column indices
        column_indices = findall(x -> x, valid_columns[:])

        # Check if any columns remain
        if isempty(column_indices)
            println("Empty community detected, returning a metric with NaN values.")
            DNCI_result = DataFrame(group1=group_combinations[i][1], group2=group_combinations[i][2], 
                        DNCI=NaN, CI_DNCI=NaN, S_DNCI=NaN,
                        status="empty_community") #Case 2
            return DNCI_result
        end
    
        # Subset the matrix using these indices
        paired_x = paired_x[:, column_indices]

        # Ensure that the paired_x matrix has at least two species
        if size(paired_x, 2) == 1
            println("Only one species present, not enough species in the selected groups to calculate DNCI. Skipping this group pair.")
            DNCI_result = DataFrame(group1=group_combinations[i][1], group2=group_combinations[i][2], 
                          DNCI=NaN, CI_DNCI=NaN, S_DNCI=NaN,
                          status= "only_one_species_exists") #Case 3
            return DNCI_result
        end


        group_pair = vcat(
            fill(group_combinations[i][1], size(splitx[group_combinations[i][1]], 1)),
            fill(group_combinations[i][2], size(splitx[group_combinations[i][2]], 1)))

        #Calculate DNCI for the group pair using the internal function
        DNCI_result = Internal.DNCI_ses(paired_x, group_pair, Nperm; count=Nperm_count)
            
        append!(ddelta, DNCI_result, promote=true)

        GC.gc()
    end
    return ddelta
end
