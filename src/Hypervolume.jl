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
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0              -0.332365                 -0.189471
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])
143×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │              -1.42803                  2.50661
   2 │               0.749608                 0.105576
   3 │               0.526123                 0.390608
   4 │              -1.20952                  2.1533
   5 │               0.855213                 0.619396
  ⋮  │           ⋮                        ⋮
 139 │              -0.466708                -0.207561
 140 │              -1.13481                 -0.447443
 141 │               1.11319                 -0.779452
 142 │               0.139825                 1.46642
 143 │              -0.3608                  -1.10642
                                        133 rows omitted  

julia> result = MVNH_det(data; var_names=["Temperature", "Precipitation"])
1×4 DataFrame
 Row │ total    correlation  Temperature  Precipitation 
     │ Float64  Float64      Float64      Float64       
─────┼──────────────────────────────────────────────────
   1 │ 1.01755     0.975263      1.03353         1.0095
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
    df_data = DataFrame(total = det(COV), correlation = rho)

    # Add each variance with corresponding variable name
    for (i, var) in enumerate(var_names)
        df_data[!, var] .= s[i]
    end

    # Convert to a DataFrame
    result = df_data

    return result
end

"""
    MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) -> DataFrame

Calculate niche dissimilarity between two species based on their environmental variables, using the Bhattacharyya distance and its components.

Arguments
- `data_1::DataFrame`: DataFrame for the first species, where each row represents an observation (presence only, need to filter out absences) and columns represent environmental variables.
- `data_2::DataFrame`: DataFrame for the second species, with the same structure as `data_1`.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names "variable1", "variable2", etc. will be used.

Returns
- `DataFrame`: A dataframe containing three metrics and their components:
  - `"Bhattacharyya_distance"`: The total Bhattacharyya distance and its components
  - `"Mahalanobis_distance"`: The Mahalanobis component of the Bhattacharyya distance
  - `"Determinant_ratio"`: The determinant ratio component of the Bhattacharyya distance

  Each metric contains:
  - `total`: The total value of the respective distance measure
  - `correlation`: The correlation component of the distance measure
  - One value for each environmental variable showing its contribution to the distance

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
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0              -0.332365                 -0.189471
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48725 rows omitted


julia> data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

143×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │              -1.42803                  2.50661
   2 │               0.749608                 0.105576
   3 │               0.526123                 0.390608
   4 │              -1.20952                  2.1533
   5 │               0.855213                 0.619396
  ⋮  │           ⋮                        ⋮
 139 │              -0.466708                -0.207561
 140 │              -1.13481                 -0.447443
 141 │               1.11319                 -0.779452
 142 │               0.139825                 1.46642
 143 │              -0.3608                  -1.10642
                                        133 rows omitted

julia> data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])
58×2 DataFrame
 Row │ normalized_temperature  normalized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │              0.72785                  0.193376
   2 │             -1.66283                 -0.135811
   3 │             -0.384102                 1.78937
   4 │             -1.74131                  0.695404
   5 │              0.847836                 0.718538
  ⋮  │           ⋮                        ⋮
  54 │             -2.42433                  0.492283
  55 │              1.2777                   1.69026
  56 │              0.128285                 0.468629
  57 │             -1.39241                  0.762875
  58 │              0.617823                -0.50501
                                         48 rows omitted
                   
julia> result = MVNH_dissimilarity(data_1, data_2; var_names=["Temperature", "Precipitation"])
3×5 DataFrame
 Row │ metric                  total      correlation  Temperature  Precipitation 
     │ String                  Float64    Float64      Float64      Float64       
─────┼────────────────────────────────────────────────────────────────────────────
   1 │ Bhattacharyya_distance  0.0213524   0.0073495    0.00345875     0.0105442
   2 │ Mahalanobis_distance    0.0106699  -0.00037229   0.00342011     0.00762205
   3 │ Determinant_ratio       0.0106826   0.00772179   3.86384e-5     0.00292214
```
"""
function MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(data_1, 2)]
    end

    # Compute covariance matrices for both datasets
    db1 = Matrix(data_1)
    db2 = Matrix(data_2)

    # Get dimensions
    n1, p = size(db1)
    n2 = size(db2, 1)

    # Compute the means for both datasets and convert them to vectors
    u1 = vec(mean(db1, dims=1)) # Convert 1-row matrix to vector
    u2 = vec(mean(db2, dims=1))
    # Compute covariance matrices
    S1 = ((db1 .- u1')' * (db1 .- u1')) / (n1 - 1)
    S2 = ((db2 .- u2')' * (db2 .- u2')) / (n2 - 1)

    # Compute average covariance
    S_avg = (S1 + S2) / 2 

    # Extract diagonal elements once for reuse
    diag_S1 = diag(S1)
    diag_S2 = diag(S2)
    diag_S_avg = diag(S_avg)

    # Calculate the Mahalanobis component
    diff = u1 .- u2  # Element-wise subtraction of vectors


    # Define mahalanobis_dist before try/catch
    mahalanobis_dist = 0.0
    # Use Cholesky decomposition 
        try
            C = cholesky(Symmetric(S_avg))
            v = C.L \ diff
            mahalanobis_dist = (1/8) * dot(v, v)
        catch
            # Fallback to direct method if Cholesky fails
            mahalanobis_dist = (1/8) * (diff' * (S_avg \ diff))
        end

    # Calculate 1D components
    mahalanobis_1d = (1 / 8) * (diff.^2) ./ diag(S_avg)
    
    # Calculate the correlation component
    mahalanobis_cor = mahalanobis_dist - sum(mahalanobis_1d)

    # Calculate the determinant ratio component
    logdet_S1 = logdet(Symmetric(S1))
    logdet_S2 = logdet(Symmetric(S2))
    logdet_S_avg = logdet(Symmetric(S_avg))

    determinant_ratio = (1/2) * (logdet_S_avg - 0.5 * (logdet_S1 + logdet_S2))
    
    # Calculate 1D components for determinant
    determinant_1d = (1/2) * log.(diag_S_avg ./ sqrt.(diag_S1 .* diag_S2))
    determinant_cor = determinant_ratio - sum(determinant_1d)

    # Total Bhattacharyya distance
    bhattacharyya_dist = mahalanobis_dist + determinant_ratio
    bhattacharyya_1d = mahalanobis_1d + determinant_1d
    bhattacharyya_cor = mahalanobis_cor + determinant_cor

    # Create a DataFrame with components as columns
    result_df = DataFrame()
    
    # Set metrics as the row names (first column)
    result_df.metric = ["Bhattacharyya_distance", "Mahalanobis_distance", "Determinant_ratio"]
       
    # Add a column for each component
    result_df.total = [bhattacharyya_dist, mahalanobis_dist, determinant_ratio]
    result_df.correlation = [bhattacharyya_cor, mahalanobis_cor, determinant_cor]
       
    # Add columns for each variable
    for (i, var) in enumerate(var_names)
        result_df[!, var] = [bhattacharyya_1d[i], mahalanobis_1d[i], determinant_1d[i]]
    end
       
    return result_df

end

"""
    average_MVNH_det(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; 
                     var_names::Vector{String}=String[]) -> Float64

Calculate the average niche hypervolume across multiple species in a community dataset.

Arguments
- `data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::Vector{String}`: Vector containing species identifiers corresponding to each observation in `data`, which must be a vector of strings.
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
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0              -0.332365                 -0.189471
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
           select(_, [:normalized_temperature, :normalized_precipitation])
           
48735×2 DataFrame
   Row │ normalized_temperature  normalized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │              0.838777                 -0.290705
     2 │             -1.10913                  -0.959396
     3 │              0.313343                 -0.660172
     4 │              0.255048                 -0.821056
     5 │             -0.402463                 -0.925731
   ⋮   │           ⋮                        ⋮
 48731 │             -0.332365                 -0.189471
 48732 │              0.516463                 -0.887027
 48733 │              0.617823                 -0.50501
 48734 │              0.391502                 -0.834642
 48735 │              0.172865                 -0.280639

julia> result = average_MVNH_det(data, df.Presence, String.(df.Species); var_names=["Temperature", "Precipitation"])
1.0026800359830965
```                                                                                                                     
"""
function average_MVNH_det(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; var_names::Vector{String}=String[])
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
            filter(row -> row.species == sp && row.presence_absence > 0, _) |>
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
    average_MVNH_dissimilarity(data::DataFrame, presence_absence::Vector{Int}, species::Union{Integer, String, PooledArrays.PooledVector}; 
                              var_names::Vector{String}=String[]) -> Float64

Calculate the average niche dissimilarity between all unique pairs of species in a community dataset using Bhattacharyya distance.

Arguments
- `data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::Vector{String}`: Vector containing species identifiers corresponding to each observation in `data`, which must be a vector of strings.
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
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.838777                 -0.290705
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.10913                  -0.959396
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               0.313343                 -0.660172
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0               0.255048                 -0.821056
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0              -0.402463                 -0.925731
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0              -0.332365                 -0.189471
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.516463                 -0.887027
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               0.617823                 -0.50501
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.391502                 -0.834642
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.172865                 -0.280639
                                                                                                                                            48725 rows omitted

julia> data = @pipe df |> 
           select(_, [:normalized_temperature, :normalized_precipitation])    
48735×2 DataFrame
   Row │ normalized_temperature  normalized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │              0.838777                 -0.290705
     2 │             -1.10913                  -0.959396
     3 │              0.313343                 -0.660172
     4 │              0.255048                 -0.821056
     5 │             -0.402463                 -0.925731
   ⋮   │           ⋮                        ⋮
 48731 │             -0.332365                 -0.189471
 48732 │              0.516463                 -0.887027
 48733 │              0.617823                 -0.50501
 48734 │              0.391502                 -0.834642
 48735 │              0.172865                 -0.280639
                                        48725 rows omitted

julia> result = average_MVNH_dissimilarity(data, df.Presence, String.(df.Species); var_names=["Temperature", "Precipitation"])     
0.060650558475890445
```
"""
function average_MVNH_dissimilarity(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; var_names::Vector{String}=String[])
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
    
    # Return average (condition ? value_if_true : value_if_false)
    return pair_count > 0 ? total_dissimilarity / pair_count : 0.0
end