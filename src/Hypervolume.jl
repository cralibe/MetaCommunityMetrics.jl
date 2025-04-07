# src/Hypervolume.jl


"""
    MVNH_det(data::DataFrame; var_names::Vector{String}=String[]) -> DataFrame

Calculate the niche hypervolume of a species based on environmental variables.

Arguments
- `data::DataFrame`: DataFrame where each row represents an observation of a species (presence only, need to filter out absences) and columns represent environmental variables.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names "variable1", "variable2", etc. will be used.

Returns
- `DataFrame`: A DataFrame containing:
  - `Correlation`: The correlation component (calculated as det(COV)/prod(variances))
  - One column for each environmental variable showing its variance
  - `total`: The total hypervolume (calculated as the determinant of the covariance matrix)

Details
- Environmental variables are assumed to follow a multivariate normal distribution, otherwise transformation to normal distribution is recommended before using this function.
- Variables should be normalized before using this function to avoid bias from different scales
- The function computes the covariance matrix of the input data, extracts variances, and calculates the determinant
- This function is a Julia implementation of the `MVNH_det` function from the R package `MVNH` (GPL-3)
- Original package and documentation: https://github.com/lvmuyang/MVNH

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics, UnicodePlots

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.80987                  -0.290381
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.12523                   0.750317
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5              -1.10775                  -1.87583
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0              -0.418417                  0.0964911
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0               0.287892                 -0.0272079
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               0.143276                  2.37981
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.148338                  1.4683
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               1.01169                  -0.485298
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.0284359                -0.392446
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.170814                  1.14892
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])
143×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │             -1.40523                   0.746757
   2 │             -2.61504                   0.229756
   3 │              0.199089                  0.347726
   4 │             -3.02404                  -0.166723
   5 │              0.835775                  0.556397
  ⋮  │           ⋮                        ⋮
 139 │              0.573287                  0.79109
 140 │             -0.0933359                -2.29843
 141 │              0.510836                 -1.97484
 142 │             -0.954948                 -0.379064
 143 │             -0.194825                 -0.167562
                                        133 rows omitted    

julia> UnicodePlots.histogram(data.normalized_temperature, 
                         nbins = 10, 
                         title = "Temperature Histogram",  
                         width=20,  
                         height=5)   
                 Temperature Histogram 
                ┌                    ┐ 
   [-4.0, -3.0) ┤▍ 1                   
   [-3.0, -2.0) ┤█▎ 4                  
   [-2.0, -1.0) ┤████▉ 16              
   [-1.0,  0.0) ┤███████████████▏ 49   
   [ 0.0,  1.0) ┤████████████████  52  
   [ 1.0,  2.0) ┤█████▊ 19             
   [ 2.0,  3.0) ┤▋ 2                   
                └                    ┘ 
                       Frequency  

julia> UnicodePlots.histogram(data.normalized_precipitation, 
                                nbins = 10, 
                                title = "Precipitation Histogram",  
                                width=20,  
                                height=5)  
                Precipitation Histogram 
                ┌                    ┐ 
   [-2.5, -2.0) ┤██  4                 
   [-2.0, -1.5) ┤██▌ 5                 
   [-1.5, -1.0) ┤█████  10             
   [-1.0, -0.5) ┤████████▌ 17          
   [-0.5,  0.0) ┤████████████████  32  
   [ 0.0,  0.5) ┤████████████████  32  
   [ 0.5,  1.0) ┤████████████  24      
   [ 1.0,  1.5) ┤██▌ 5                 
   [ 1.5,  2.0) ┤████▌ 9               
   [ 2.0,  2.5) ┤██▌ 5                 
                └                    ┘ 
                       Frequency   

julia> result = MVNH_det(data; var_names=["Temperature", "Precipitation"])
1×4 DataFrame
 Row │ Correlation  Precipitation  Temperature  total    
     │ Float64      Float64        Float64      Float64  
─────┼───────────────────────────────────────────────────
   1 │    0.999758       0.942899      0.99626  0.939145

```
"""
function MVNH_det(data::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(data, 2)]
    end

    # Convert DataFrame to Matrix (drop non-numeric columns if any)
    data_matrix = Matrix(data[:, :])

    # Compute the covariance matrix from the data
    COV = cov(data_matrix)

    # Calculate the variance (diagonal of the covariance matrix)
    s = diag(COV)

    # Calculate the correlation component (determinant ratio)
    rho = det(COV) / prod(s)

    # Prepare the data for creating the DataFrame
    # Start with the determinant and correlation component
    df_data = Dict("total" => det(COV), "Correlation" => rho)

    # Add each variance with corresponding variable name
    for i in 1:length(var_names)
        df_data[var_names[i]] = s[i]
    end

    # Convert to a DataFrame
    result = DataFrame(df_data)

    return result
end

"""
    MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) -> Dict{String, DataFrame}

Calculate niche dissimilarity between two species based on their environmental variables, using the Bhattacharyya distance and its components.

Arguments
- `data_1::DataFrame`: DataFrame for the first species, where each row represents an observation (presence only, need to filter out absences) and columns represent environmental variables.
- `data_2::DataFrame`: DataFrame for the second species, with the same structure as `data_1`.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names "variable1", "variable2", etc. will be used.

Returns
- `Dict{String, DataFrame}`: A dictionary containing three DataFrames:
  - `"Bhattacharyya_distance"`: The total Bhattacharyya distance and its components
  - `"Mahalanobis_distance"`: The Mahalanobis component of the Bhattacharyya distance
  - `"Determinant_ratio"`: The determinant ratio component of the Bhattacharyya distance

  Each DataFrame contains:
  - `total`: The total value of the respective distance measure
  - `correlation`: The correlation component of the distance measure
  - One column for each environmental variable showing its contribution to the distance

# Details
- The Bhattacharyya distance is calculated as the sum of two components:
  1. Mahalanobis component: (1/8) × (μ₁-μ₂)ᵀ × (S₁+S₂)/2⁻¹ × (μ₁-μ₂)
  2. Determinant ratio component: (1/2) × log(det((S₁+S₂)/2) / sqrt(det(S₁) × det(S₂)))
- Each component is further decomposed into individual variable contributions and correlation effects
- Environmental variables are assumed to follow a multivariate normal distribution, otherwise transformation to normal distribution is recommended before using this function.
- Variables should be normalized before using this function to avoid bias from different scales
- This function is a Julia implementation inspired by the `MVNH` R package (GPL-3)
- Original package and documentation: https://github.com/lvmuyang/MVNH

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.80987                  -0.290381
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.12523                   0.750317
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5              -1.10775                  -1.87583
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0              -0.418417                  0.0964911
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0               0.287892                 -0.0272079
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               0.143276                  2.37981
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.148338                  1.4683
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               1.01169                  -0.485298
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.0284359                -0.392446
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.170814                  1.14892
                                                                                                                                            48725 rows omitted


julia> data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])
143×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │             -1.40523                   0.746757
   2 │             -2.61504                   0.229756
   3 │              0.199089                  0.347726
   4 │             -3.02404                  -0.166723
   5 │              0.835775                  0.556397
  ⋮  │           ⋮                        ⋮
 139 │              0.573287                  0.79109
 140 │             -0.0933359                -2.29843
 141 │              0.510836                 -1.97484
 142 │             -0.954948                 -0.379064
 143 │             -0.194825                 -0.167562
                                        133 rows omitted

julia> data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])
58×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │               0.319133                -1.26372
   2 │               0.227851                -1.30392
   3 │               0.688124                 0.63633
   4 │              -0.630506                 0.269884
   5 │              -0.600362                -0.601565
  ⋮  │           ⋮                        ⋮
  54 │              -0.560579                -1.1251
  55 │              -0.216268                 0.673043
  56 │               1.04894                 -0.264528
  57 │              -1.52279                  1.80578
  58 │               1.01169                 -0.485298
                                         48 rows omitted
                   
julia> result = MVNH_dissimilarity(data_1, data_2; var_names=["Temperature", "Precipitation"])
Dict{String, DataFrame} with 3 entries:
  "Determinant_ratio"      => 1×4 DataFrame…
  "Bhattacharyya_distance" => 1×4 DataFrame…
  "Mahalanobis_distance"   => 1×4 DataFrame…

julia> result["Determinant_ratio"]
1×4 DataFrame
 Row │ total      correlation  Temperature  Precipitation 
     │ Float64    Float64      Float64      Float64       
─────┼────────────────────────────────────────────────────
   1 │ 0.0048021  0.000767901  0.000672539     0.00336166

julia> result["Bhattacharyya_distance"]
1×4 DataFrame
 Row │ total       correlation  Temperature  Precipitation 
     │ Float64     Float64      Float64      Float64       
─────┼─────────────────────────────────────────────────────
   1 │ 0.00932099   0.00102573   0.00296459     0.00533067

julia> result["Mahalanobis_distance"]
1×4 DataFrame
 Row │ total       correlation  Temperature  Precipitation 
     │ Float64     Float64      Float64      Float64       
─────┼─────────────────────────────────────────────────────
   1 │ 0.00451889  0.000257828   0.00229205     0.00196901

```
"""
function MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(db1, 2)]
    end

    # Compute covariance matrices for both datasets
    db1 = Matrix(data_1[:, :])
    db2 = Matrix(data_2[:, :])
    S1 = cov(db1)
    S2 = cov(db2)
    S_avg = (S1 + S2) / 2  # Average covariance matrix

    # Compute the means for both datasets and convert them to vectors
    u1 = vec(mean(db1, dims=1))  # Convert 1-row matrix to vector
    u2 = vec(mean(db2, dims=1))  # Convert 1-row matrix to vector

    # Calculate the Mahalanobis component
    diff = u1 .- u2  # Element-wise subtraction of vectors
    mahalanobis_dist = (1 / 8) * (diff' * inv(S_avg) * diff)
    mahalanobis_1d = (1 / 8) * (diff.^2) ./ diag(S_avg)
    mahalanobis_cor = mahalanobis_dist[1] - sum(mahalanobis_1d)

    # Calculate the determinant ratio component
    determinant_ratio = (1 / 2) * log(det(S_avg) / sqrt(det(S1) * det(S2)))
    determinant_1d = (1 / 2) * log.(diag(S_avg) ./ sqrt.(diag(S1) .* diag(S2)))
    determinant_cor = determinant_ratio - sum(determinant_1d)

    # Total Bhattacharyya distance
    bhattacharyya_dist = mahalanobis_dist[1] + determinant_ratio
    bhattacharyya_1d = mahalanobis_1d + determinant_1d
    bhattacharyya_cor = mahalanobis_cor + determinant_cor

    # Create Dictionary for each component
    mahalanobis_results = Dict("total" => mahalanobis_dist[1], "correlation" => mahalanobis_cor)
    determinant_results = Dict("total" => determinant_ratio, "correlation" => determinant_cor)
    bhattacharyya_results = Dict("total" => bhattacharyya_dist, "correlation" => bhattacharyya_cor)

    for i in 1:length(var_names)
        mahalanobis_results[var_names[i]] = mahalanobis_1d[i]
        determinant_results[var_names[i]] = determinant_1d[i]
        bhattacharyya_results[var_names[i]] = bhattacharyya_1d[i]
    end

    # Convert to DataFrames
    mahalanobis_df = DataFrame(mahalanobis_results)
    determinant_df = DataFrame(determinant_results)
    bhattacharyya_df = DataFrame(bhattacharyya_results)

    # Set the desired column order: ["total", "correlation", var_names...]
    column_order = vcat(["total", "correlation"], var_names)

    # Reorder columns in the DataFrames
    mahalanobis_df = select(mahalanobis_df, column_order)
    determinant_df = select(determinant_df, column_order)
    bhattacharyya_df = select(bhattacharyya_df, column_order)

    # Return all three DataFrames in a dictionary
    return Dict(
        "Bhattacharyya_distance" => bhattacharyya_df,
        "Mahalanobis_distance" => mahalanobis_df,
        "Determinant_ratio" => determinant_df
    )

end

"""
    average_MVNH_det(data::DataFrame, presence_absence::Vector{Int}, species::Union{AbstractVector, String}; 
                     var_names::Vector{String}=String[]) -> Float64

Calculate the average niche hypervolume across multiple species in a community dataset.

Arguments
- `data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::Union{AbstractVector, String}`: Vector containing species identifiers corresponding to each observation in `data`.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names will be used.

Returns
- `Float64`: The average hypervolume across all species with presence data.

Details
- For each unique species, the function:
  1. Filters observations where the species is present (presence_absence > 0)
  2. Calculates the niche hypervolume using the `MVNH_det` function
  3. Extracts the total hypervolume value
- The function then computes the mean of all individual species hypervolumes
- Species with no presence data are skipped in the calculation
- Environmental variables are assumed to follow a multivariate normal distribution, otherwise transformation to normal distribution is recommended before using this function.
- Variables should be normalized before using this function to avoid bias from different scales

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.80987                  -0.290381
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.12523                   0.750317
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5              -1.10775                  -1.87583
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0              -0.418417                  0.0964911
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0               0.287892                 -0.0272079
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               0.143276                  2.37981
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.148338                  1.4683
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               1.01169                  -0.485298
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.0284359                -0.392446
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.170814                  1.14892
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
           select(_, [:normalized_temperature, :normalized_precipitation])
           
48735×2 DataFrame
   Row │ normalized_temperature  normalized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │              0.80987                  -0.290381
     2 │             -1.12523                   0.750317
     3 │             -1.10775                  -1.87583
     4 │             -0.418417                  0.0964911
     5 │              0.287892                 -0.0272079
   ⋮   │           ⋮                        ⋮
 48731 │              0.143276                  2.37981
 48732 │              0.148338                  1.4683
 48733 │              1.01169                  -0.485298
 48734 │              0.0284359                -0.392446
 48735 │              0.170814                  1.14892
                                        48725 rows omitted

julia> result = average_MVNH_det(data, df.Presence, df.Species; var_names=["Temperature", "Precipitation"])
0.9842468737598974

```                                                                                                                     
"""
function average_MVNH_det(data::DataFrame, presence_absence::Vector{Int}, species::Union{AbstractVector, String}; var_names::Vector{String}=String[])
    # Ensure the presence_absence and species vectors match the DataFrame size
    @assert length(presence_absence) == nrow(data) "Mismatch in data length"
    @assert length(species) == nrow(data) "Mismatch in species length"

    # Add presence_absence and species as columns to the DataFrame
    data = hcat(data, DataFrame(presence_absence=presence_absence, species=species))

    # Initialize final result DataFrame
    final_result = DataFrame(species=String[], hypervolume=Float64[])

    # Loop over each unique species
    for sp in unique(species)
        # Filter rows for the current species with presence > 0
        filtered_data = @pipe data |>
            filter(row -> row.species == sp && row.presence_absence == 1, _) |>
            select(_, Not(:species,:presence_absence))

        # Skip if no valid rows for the species
        if nrow(filtered_data) == 0
            continue
        end     
            
        # Compute the MVNH_det result
        result = MVNH_det(filtered_data; var_names)

        # Add the species and its total hypervolume to the final result DataFrame
        push!(final_result, (species=sp, hypervolume=result.total[1]))

    end

    # Compute the average hypervolume across species
    averaged_hypervolume = mean(final_result.hypervolume)

    return averaged_hypervolume

end

"""
    average_MVNH_dissimilarity(data::DataFrame, presence_absence::Vector{Int}, species::Union{AbstractVector, String}; 
                              var_names::Vector{String}=String[]) -> Float64

Calculate the average niche dissimilarity between all unique pairs of species in a community dataset using Bhattacharyya distance.

Arguments
- `data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::Union{AbstractVector, String}`: Vector containing species identifiers corresponding to each observation in `data`.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names will be used.

Returns
- `Float64`: The average Bhattacharyya distance across all unique species pairs.

Details
- For each unique pair of species, the function:
  1. Filters observations where each species is present (presence_absence > 0)
  2. Calculates the niche dissimilarity using the `MVNH_dissimilarity` function
  3. Extracts the total Bhattacharyya distance value
- The function then computes the mean of all pairwise Bhattacharyya distances
- Species pairs where either species has no presence data are skipped
- Each species pair is processed only once (i.e., sp1-sp2 is calculated, but sp2-sp1 is skipped)
- Environmental variables are assumed to follow a multivariate normal distribution
- Variables should be normalized before using this function to avoid bias from different scales

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.80987                  -0.290381
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.12523                   0.750317
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5              -1.10775                  -1.87583
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0              -0.418417                  0.0964911
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0               0.287892                 -0.0272079
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               0.143276                  2.37981
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.148338                  1.4683
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               1.01169                  -0.485298
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.0284359                -0.392446
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.170814                  1.14892
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
           select(_, [:normalized_temperature, :normalized_precipitation])    
48735×2 DataFrame
   Row │ normalized_temperature  normalized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │              0.80987                  -0.290381
     2 │             -1.12523                   0.750317
     3 │             -1.10775                  -1.87583
     4 │             -0.418417                  0.0964911
     5 │              0.287892                 -0.0272079
   ⋮   │           ⋮                        ⋮
 48731 │              0.143276                  2.37981
 48732 │              0.148338                  1.4683
 48733 │              1.01169                  -0.485298
 48734 │              0.0284359                -0.392446
 48735 │              0.170814                  1.14892
                                        48725 rows omitted

julia> result = average_MVNH_dissimilarity(data, df.Presence, df.Species; var_names=["Temperature", "Precipitation"])     
0.02923266035138391

```
"""
function average_MVNH_dissimilarity(data::DataFrame, presence_absence::Vector{Int}, species::Union{AbstractVector, String}; var_names::Vector{String}=String[])
    @assert length(presence_absence) == nrow(data) "Mismatch in data length"
    @assert length(species) == nrow(data) "Mismatch in species length"

    # Add presence_absence and species as columns to the DataFrame
    data = hcat(data, DataFrame(presence_absence=presence_absence, species=species))

    # Initialize final result DataFrame
    final_result = DataFrame(species_1=String[], species_2=String[], hypervolume_dissimilarity=Float64[])

     # Loop over each unique species pair
    for sp1 in unique(species)
        for sp2 in unique(species)
            # Skip if the same species or if the pair has already been processed or if the pair is the reverse
            if sp1 == sp2 || nrow(filter(row -> row.species_1 == sp1 && row.species_2 == sp2, final_result)) > 0 || nrow(filter(row -> row.species_1 == sp2 && row.species_2 == sp1, final_result)) > 0
                continue
            end

            # Filter rows for the current species pair with presence > 0
            filtered_data_1 = @pipe data |>
                filter(row -> row.species == sp1 && row.presence_absence > 0, _) |>
                select(_, Not(:species,:presence_absence))

            filtered_data_2 = @pipe data |>
                filter(row -> row.species == sp2 && row.presence_absence > 0, _) |>
                select(_, Not(:species,:presence_absence))

            # Skip if no valid rows for the species pair
            if nrow(filtered_data_1) == 0 || nrow(filtered_data_2) == 0
                continue
            end

            # Compute the MVNH_dissimilarity result
            result = MVNH_dissimilarity(filtered_data_1, filtered_data_2; var_names)

            # Add the species pair and its hypervolume dissimilarity to the final result DataFrame
            push!(final_result, (species_1=sp1, species_2=sp2, hypervolume_dissimilarity=result["Bhattacharyya_distance"].total[1]))
        
        end
    end

    mean_hypervolume_dissimilarity = mean(final_result.hypervolume_dissimilarity)
            
    return mean_hypervolume_dissimilarity
    
end