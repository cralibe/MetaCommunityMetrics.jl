# src/Hypervolume.jl


"""
    MVNH_det(env_data::DataFrame; var_names::Vector{String}=String[]) -> DataFrame

Calculate the niche hypervolume of a species based on environmental variables.

Arguments
- `env_data::DataFrame`: DataFrame where each row represents an observation of a species (presence only, need to filter out absences) and columns represent environmental variables.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names "variable1", "variable2", etc. will be used.

Returns
`DataFrame`: A DataFrame containing:
    - `Correlation`: The correlation component (calculated as det(COV)/prod(variances))
    - One column for each environmental variable showing its variance
    - `total`: The total hypervolume (calculated as the determinant of the covariance matrix)

Details
- Environmental variables are assumed to follow a multivariate normal distribution, otherwise transformation to normal distribution is recommended before using this function.
- Variables should be standardized before using this function to avoid bias from different scales
- The function computes the covariance matrix of the input data, extracts variances, and calculates the determinant
- This function is a Julia implementation of the `MVNH_det()` function from the R package `MVNH` (GPL-3)
- Original package and documentation: https://github.com/lvmuyang/MVNH

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics, UnicodePlots

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

julia> env_data = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:standardized_temperature, :standardized_precipitation])
143×2 DataFrame
 Row │ standardized_temperature  standardized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │            -0.37813                    0.13009
   2 │            -0.00856861                 0.237183
   3 │            -0.664638                  -0.772406
   4 │             2.05431                    0.451875
   5 │            -0.39968                   -0.719024
  ⋮  │           ⋮                        ⋮
 139 │             1.85574                   -0.583737
 140 │             0.0953878                  1.21099
 141 │            -1.02227                    1.33501
 142 │            -0.400246                  -0.438892
 143 │            -0.817817                   0.418038
                                        133 rows omitted

julia> result = MVNH_det(env_data; var_names=["Temperature", "Precipitation"])
1×4 DataFrame
 Row │ total    correlation  Temperature  Precipitation 
     │ Float64  Float64      Float64      Float64       
─────┼──────────────────────────────────────────────────
   1 │ 1.15268     0.999732     0.962495        1.19792
```
"""
function MVNH_det(env_data::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(env_data, 2)]
    end

    # Convert DataFrame to Matrix (drop non-numeric columns if any)
    data_matrix = Matrix(env_data[:, :])

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
    MVNH_dissimilarity(env_data_1::DataFrame, env_data_2::DataFrame; var_names::Vector{String}=String[]) -> DataFrame

Calculate niche dissimilarity between two species based on their environmental variables, using the Bhattacharyya distance and its components.

Arguments
- `env_data_1::DataFrame`: DataFrame for the first species, where each row represents an observation (presence only, need to filter out absences) and columns represent environmental variables.
- `env_data_2::DataFrame`: DataFrame for the second species, with the same structure as `data_1`.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names "variable1", "variable2", etc. will be used.

Returns
`DataFrame`: A dataframe containing three metrics and their components:
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
- Variables should be standardized_ before using this function to avoid bias from different scales
- This function is a Julia implementation inspired by the `MVNH` R package (GPL-3)
- Original package and documentation: https://github.com/lvmuyang/MVNH

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

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


julia> env_data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:standardized_temperature, :standardized_precipitation])

143×2 DataFrame
 Row │ standardized_temperature  standardized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │            -0.37813                    0.13009
   2 │            -0.00856861                 0.237183
   3 │            -0.664638                  -0.772406
   4 │             2.05431                    0.451875
   5 │            -0.39968                   -0.719024
  ⋮  │           ⋮                        ⋮
 139 │             1.85574                   -0.583737
 140 │             0.0953878                  1.21099
 141 │            -1.02227                    1.33501
 142 │            -0.400246                  -0.438892
 143 │            -0.817817                   0.418038
                                        133 rows omitted

julia> env_data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:standardized_temperature, :standardized_precipitation])
58×2 DataFrame
 Row │ standardized_temperature  standardized_precipitation 
     │ Float64                 Float64                  
─────┼──────────────────────────────────────────────────
   1 │              -0.229864                  1.84371
   2 │               0.460218                 -0.624328
   3 │              -1.03283                  -1.16451
   4 │               0.675006                 -0.120586
   5 │               0.40729                   0.20034
  ⋮  │           ⋮                        ⋮
  54 │              -0.870299                 -0.235392
  55 │               0.504555                 -1.50887
  56 │               2.03065                  -0.740789
  57 │              -0.174396                  0.448461
  58 │               0.547169                  1.03257
                                         48 rows omitted
                   
julia> result = MVNH_dissimilarity(env_data_1, env_data_2; var_names=["Temperature", "Precipitation"])
3×5 DataFrame
 Row │ metric                  total       correlation  Temperature  Precipitation 
     │ String                  Float64     Float64      Float64      Float64       
─────┼─────────────────────────────────────────────────────────────────────────────
   1 │ Bhattacharyya_distance  0.00980771  0.00015205    0.00388058     0.00577508
   2 │ Mahalanobis_distance    0.00664862  5.06232e-6    0.00234902     0.00429454
   3 │ Determinant_ratio       0.00315908  0.000146988   0.00153156     0.00148054
```
"""
function MVNH_dissimilarity(env_data_1::DataFrame, env_data_2::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(env_data_1, 2)]
    end

    # Compute covariance matrices for both datasets
    db1 = Matrix(env_data_1)
    db2 = Matrix(env_data_2)

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
    average_MVNH_det(env_data::DataFrame, presence_absence::Vector{Int}, species::AbstractVector; 
                     var_names::Vector{String}=String[]) -> Float64

Calculate the average niche hypervolume across multiple species in a community dataset.

Arguments
- `env_data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::AbstractVector`: Vector containing species identifiers corresponding to each observation in `data`, which must be a vector of strings.
- `var_names::Vector{String}=String[]`: Optional vector specifying names for the environmental variables. If empty, default names will be used.

Returns
- `Float64`: The average hypervolume across all species with presence data.

Details
- For each unique species, the function:
  1. Filters observations where the species is present (presence_absence > 0)
  2. Calculates the niche hypervolume using the `MVNH_det` function
  3. Extracts the total hypervolume value
- The function then computes the mean of all individual species hypervolumes.
- Species with no presence data are skipped in the calculation.
- Environmental variables are assumed to follow a multivariate normal distribution, otherwise transformation to normal distribution is recommended before using this function.
- User should standardize the environmental variables before using this function to avoid bias from different scales. So user should transform them to normal distribution when necessary.
- User should remove singletons (species only occupied one site) before using this function to avoid NaN values in the hypervolume calculation.

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = @pipe load_sample_data() |>
                      groupby(_, :Species) |>
                      filter(row -> sum(row.Presence) > 1, _)|>
                      DataFrame(_)
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

julia> env_data = @pipe df |> 
           select(_, [:standardized_temperature, :standardized_precipitation]) 
53352×2 DataFrame
   Row │ standardized_temperature  standardized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │               0.829467              -1.4024
     2 │              -1.12294               -0.0519895
     3 │              -0.409808              -0.803663
     4 │              -1.35913               -0.646369
     5 │               0.0822                 1.09485
   ⋮   │           ⋮                        ⋮
 53348 │              -0.571565              -0.836345
 53349 │              -2.33729               -0.398522
 53350 │               0.547169               1.03257
 53351 │              -0.815015               0.95971
 53352 │               0.48949               -1.59416
                                        53342 rows omitted

julia> result = average_MVNH_det(env_ata, df.Presence, df.Species; var_names=["Temperature", "Precipitation"])
1.2103765096417536
```                                                                                                                     
"""
function average_MVNH_det(env_data::DataFrame, presence_absence::Vector{Int}, species::AbstractVector; var_names::Vector{String}=String[])
    # Ensure the presence_absence and species vectors match the DataFrame size
    @assert length(presence_absence) == nrow(env_data) "Mismatch in data length"
    @assert length(species) == nrow(env_data) "Mismatch in species length"

    # Add presence_absence and species as columns to the DataFrame
    data = hcat(env_data, DataFrame(presence_absence=presence_absence, species=species))

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
    average_MVNH_dissimilarity(env_data::DataFrame, presence_absence::Vector{Int}, species::AbstractVector; 
                              var_names::Vector{String}=String[]) -> Float64

Calculate the average niche dissimilarity between all unique pairs of species in a community dataset using Bhattacharyya distance.

Arguments
- `env_data::DataFrame`: DataFrame containing environmental variables where each row represents an observation.
- `presence_absence::Vector{Int}`: Vector indicating presence (1) or absence (0) for each observation in `data`.
- `species::AbstractVector`: Vector containing species identifiers corresponding to each observation in `data`, which must be a vector of strings.
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
- Environmental variables are assumed to follow a multivariate normal distribution. So user should transform them to normal distribution when necessary.
- User should standardize the environmental variables before using this function to avoid bias from different scales.
- User should remove singletons (species only occupied one site) before using this function to avoid NaN values in the hypervolume calculation.

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = @pipe load_sample_data() |>
                      groupby(_, :Species) |>
                      filter(row -> sum(row.Presence) > 1, _)|>
                      DataFrame(_)
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

julia> env_data = @pipe df |> 
           select(_, [:standardized_temperature, :standardized_precipitation])    
53352×2 DataFrame
   Row │ standardized_temperature  standardized_precipitation 
       │ Float64                 Float64                  
───────┼──────────────────────────────────────────────────
     1 │               0.829467              -1.4024
     2 │              -1.12294               -0.0519895
     3 │              -0.409808              -0.803663
     4 │              -1.35913               -0.646369
     5 │               0.0822                 1.09485
   ⋮   │           ⋮                        ⋮
 53348 │              -0.571565              -0.836345
 53349 │              -2.33729               -0.398522
 53350 │               0.547169               1.03257
 53351 │              -0.815015               0.95971
 53352 │               0.48949               -1.59416
                                        53342 rows omitted

julia> result = average_MVNH_dissimilarity(env_data, df.Presence, df.Species; var_names=["Temperature", "Precipitation"])     
0.03059942936454443
```
"""
function average_MVNH_dissimilarity(env_data::DataFrame, presence_absence::Vector{Int}, species::AbstractVector; var_names::Vector{String}=String[])
    @assert length(presence_absence) == nrow(env_data) "Mismatch in data length"
    @assert length(species) == nrow(env_data) "Mismatch in species length"

    # Pre-compute unique species list once
    unique_species = unique(species)
    
    # Convert DataFrame to matrix once (this is much more efficient)
    data_matrix = Matrix(env_data)
    
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