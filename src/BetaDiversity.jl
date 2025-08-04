# src/BetaDiveristy.jl

"""
    beta_diversity(mat::Matrix; quant::Bool) -> DataFrame

Calculate beta diversity for a given biodiversity data. This function supports both binary (presence/absence) and quantitative (abundance) data.

Arguments
- `mat::Matrix`: A matrix where each row represents a site and each column represents a species. The elements of the matrix should represent the presence/absence or abundance of species.
- `quant::Bool`: Specifies the data type for analysis. When `false`, treats data as binary presence/absence, converting any quantitative values and applying Jaccard-based indices. When `true`, treats data as quantitative abundance data and applies Ruzicka-based indices.

Returns
- `DataFrame`: A DataFrame with the following columns:
    - `BDtotal`: Total beta diversity, which captures the overall dissimilarity between local communities.
    - `Repl`: Replacement component of diversity, which reflects species turnover between sites: how much communities differ due to having different species compositions rather than different species counts.
    - `RichDif`: Richness difference component of diversity, which captures differences in the number of species between communities.

Details
- Empty sites (empty rows in the matrix) have to be removed before calculation.
- Species that did not occupy any sites in the data (empty columns in the matrix) have to be removed before calculation.
- This function is a translation/adaptation of the beta.dov.comp function from the R package `adespatial`,licensed under GPL-3.
- Original package and documentation available at: https://cran.r-project.org/web/packages/adespatial/index.html

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String   Int64      Int64     Float64   Float64    Float64?                Float64?                 
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

julia> matrix_with_abundance = @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> 
           filter(row -> row[:Sampling_date_order] == 1, _) |> 
           unstack(_, :Species, :Abundance, fill=0) |>  
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:] 
15×5 Matrix{Union{Missing, Int64}}:
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 2  0  0  0  0
 1  0  0  1  0
 4  0  0  1  0
 1  0  0  0  0
 ⋮           
 0  0  0  2  0
 1  0  0  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0  1  0  0
 0  0  0  0  1

julia> matrix_with_presence =  @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> 
           filter(row -> row[:Sampling_date_order] == 1, _) |>
           unstack(_, :Species, :Presence, fill=0) |> 
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix      
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:]
15×5 Matrix{Int64}:
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  1  0
 1  0  0  1  0
 1  0  0  0  0
 ⋮           
 0  0  0  1  0
 1  0  0  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0  1  0  0
 0  0  0  0  1

julia> result_using_abanduce_data_1 = beta_diversity(matrix_with_abundance; quant=true)
1×3 DataFrame
 Row │ BDtotal   Repl     RichDif  
     │ Float64   Float64  Float64  
─────┼─────────────────────────────
   1 │ 0.390317   0.2678  0.122517

julia> result_using_abanduce_data_2 = beta_diversity(matrix_with_abundance; quant=false)
1×3 DataFrame
 Row │ BDtotal   Repl      RichDif   
     │ Float64   Float64   Float64   
─────┼───────────────────────────────
   1 │ 0.357143  0.284127  0.0730159

julia> result_using_binary_data = beta_diversity(matrix_with_presence; quant=false)
1×3 DataFrame
 Row │ BDtotal   Repl      RichDif   
     │ Float64   Float64   Float64   
─────┼───────────────────────────────
   1 │ 0.357143  0.284127  0.0730159
```

"""
function beta_diversity(mat::Matrix; quant::Bool)
    #Error handling
    if isempty(mat)
        throw(ArgumentError("Input matrix is empty."))
    end

    # Error if any row has a sum of 0
    if any(sum(mat, dims=2) .== 0)
        throw(ArgumentError("One or more rows have a sum of zero, indicating no data for these rows."))
    end

    # Error if any column has a sum of 0
    if any(sum(mat, dims=1) .== 0)
        throw(ArgumentError("One or more columns have a sum of zero, indicating insufficient variation."))
    end

    
    # Check if there is variation in the data matrix
    if sum((mat .- mean(mat, dims=1)).^2) == 0
        println("The data matrix has no variation in beta space")
        return DataFrames.DataFrame(BDtotal=0, Repl=0, RichDif=0)
    end
    
    if !quant #Binary data, presence-absence data
        # Calculate Podani family, Jaccard-based indices
        mat = ifelse.(mat .> 0, 1, 0)
        a = mat * transpose(mat)
        b = mat * (1 .- transpose(mat))
        c=(1 .- mat) * transpose(mat)
        min_bc=min.(b, c)
        repl=2 .* min_bc
        rich=abs.(b - c)
            
        #Jaccard-based components
        repl = repl ./ (a + b + c)
        rich = rich ./ (a + b + c)
        D = (b + c) ./ (a + b + c)

        n = size(mat,1) #Count the number of rows
        
        total_div = (sum(D)/2) / (n * (n - 1))
        repl_div = (sum(repl)/2) / (n * (n - 1))
        rich_div = (sum(rich)/2) / (n * (n - 1))

    else # Quantitative data
            
        n = size(mat, 1)
        repl = zeros(n, n)
        rich = zeros(n, n)
        D = zeros(n, n)

        # Calculate Podani family, Ruzicka-based indices
        for i in 2:n
            for j in 1:(i - 1)
                tmp = mat[i, :] - mat[j, :]
                A = sum(min.(mat[i,:], mat[j,:]))
                B = sum(tmp[tmp .> 0])
                C = -sum(tmp[tmp .< 0])
                    
                den = A + B + C 

                repl[i,j] = 2 * (min(B,C)) / den
                rich[i,j] = abs(B - C) / den
                D[i,j] = (B + C) / den
            end
        end

        # Convert matrices to distance objects
        repl = @pipe repl|>
        [i > j ? _[i, j] : 0 for i in 1:size(_, 1), j in 1:size(_, 2)]
        rich = @pipe rich|>
        [i > j ? _[i, j] : 0 for i in 1:size(_, 1), j in 1:size(_, 2)]
        D = @pipe D|>
        [i > j ? _[i, j] : 0 for i in 1:size(_, 1), j in 1:size(_, 2)]

        # Calculate diversity metrics
        repl_div = sum(repl) / (n * (n - 1))
        rich_div = sum(rich) / (n * (n - 1))
        total_div = sum(D) / (n * (n - 1))

    end

    # Return results as a DataFrame
    return DataFrames.DataFrame(
            BDtotal=total_div, 
            Repl=repl_div, 
            RichDif=rich_div)
    
end

"""
    spatial_beta_div(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector; quant::Bool) -> DataFrame

Calculate the beta diversity decompositions of a metacommunity in space based on species abundances or presence-absences using the function `beta_diversity`.

Arguments
- `abundance::AbstractVector`: Vector representing the abundance or occurence of species.
- `time::AbstractVector`: Vector representing sampling dates.
- `site::AbstractVector`: Vector representing site names or IDs.
- `species::AbstractVector`: Vector representing species names or IDs.
- `quant::Bool`: Specifies the data type for analysis. When `false`, treats data as binary presence/absence, converting any quantitative values and applying Jaccard-based indices. When `true`, treats data as quantitative abundance data and applies Ruzicka-based indices.

Returns
- `DataFrame`: A DataFrame containing the values of total beta diversity, replacement, and richness difference components in space. Columns are `spatial_BDtotal`, `spatial_Repl`, and `spatial_RichDif`.

Details
- This function uses the `beta_diversity` function to calculate beta diversity components after aggregating .
- For binary data, the function calculates Podani family, Jaccard-based indices. 
- For quantitative data, the function calculates Podani family, Ruzicka-based indices.
- The beta diversity decompositions of a metacommunity in space is a set of metrics suggested by Guzman et al. (2022) to infer metacommunity processes.

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

julia> result_using_abanduce_data_1 = spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
1×3 DataFrame
 Row │ spatial_BDtotal  spatial_Repl  spatial_RichDif 
     │ Float64          Float64       Float64         
─────┼────────────────────────────────────────────────
   1 │        0.264822      0.121882         0.142939
        
julia> result_using_abanduce_data_2 = spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ spatial_BDtotal  spatial_Repl  spatial_RichDif 
     │ Float64          Float64       Float64         
─────┼────────────────────────────────────────────────
   1 │        0.133035       0.05976        0.0732746

julia> result_using_binary_data = spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ spatial_BDtotal  spatial_Repl  spatial_RichDif 
     │ Float64          Float64       Float64         
─────┼────────────────────────────────────────────────
   1 │        0.133035       0.05976        0.0732746
```

"""
function spatial_beta_div(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector; quant::Bool)
    #Create the required data frame
    df = @pipe DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=site,
        Species=species)

    df_matrix = @pipe df |>
        groupby(_, [:Patch, :Species]) |>
        combine(_, :N => sum => :total_N) |>
        unstack(_, :Species, :total_N, fill=0) |>
        select(_, Not(:Patch, )) |>
        Matrix(_) |>
        _[:, sum(_, dims=1)[1, :] .!= 0] |> 
        _[sum(_, dims=2)[:, 1] .!= 0,:] 
    
    #Calculate beta diversity components
    components = beta_diversity(df_matrix; quant=quant)

    spatial_beta_div_summary = DataFrames.DataFrame(
        spatial_BDtotal =  components.BDtotal,
        spatial_Repl =  components.Repl,
        spatial_RichDif =  components.RichDif)

    return spatial_beta_div_summary
end

"""
    temporal_beta_div(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector;quant::Bool) -> DataFrame

Calculate the beta diversity decompositions of a metacommunity in time based on species abundances or presence-absences using the function `beta_diversity`.

Arguments
- `abundance::AbstractVector`: Vector representing the abundance or occurence of species.
- `time::AbstractVector`: Vector representing sampling dates.
- `site::AbstractVector`: Vector representing site names or IDs.
- `species::AbstractVector`: Vector representing species names or IDs.
- `quant::Bool`: Specifies the data type for analysis. When `false`, treats data as binary presence/absence, converting any quantitative values and applying Jaccard-based indices. When `true`, treats data as quantitative abundance data and applies Ruzicka-based indices.

Returns
- `DataFrame`: A DataFrame containing the values of total beta diversity, replacement, and richness difference components in time. Columns are `temporal_BDtotal`, `temporal_Repl`, and `temporal_RichDif`.

Details
- This function uses the `beta_diversity` function to calculate beta diversity components for each site.
- For binary data, the function calculates Podani family, Jaccard-based indices. 
- For quantitative data, the function calculates Podani family, Ruzicka-based indices.
- The beta diversity decompositions of a metacommunity in time is a set of metrics suggested by Guzman et al. (2022) to infer metacommunity processes.

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

julia> result_using_abanduce_data_1 = temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
1×3 DataFrame
 Row │ temporal_BDtotal  temporal_Repl  temporal_RichDif 
     │ Float64           Float64        Float64          
─────┼───────────────────────────────────────────────────
   1 │         0.311222      0.0995483          0.211674
        
julia> result_using_abanduce_data_2 = temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ temporal_BDtotal  temporal_Repl  temporal_RichDif 
     │ Float64           Float64        Float64          
─────┼───────────────────────────────────────────────────
   1 │         0.206262      0.0693664          0.136895

julia> result_using_binary_data = temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ temporal_BDtotal  temporal_Repl  temporal_RichDif 
     │ Float64           Float64        Float64          
─────┼───────────────────────────────────────────────────
   1 │         0.206262      0.0693664          0.136895
```

"""
function temporal_beta_div(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector; quant::Bool)
    #Create the required data frame
    df = @pipe DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=site,
        Species=species)

    df_matrix = @pipe df |>
        groupby(_, [:Time, :Species]) |>
        combine(_, :N => sum => :total_N) |>
        unstack(_, :Species, :total_N, fill=0) |>
        select(_, Not(:Time)) |>
        Matrix(_) |>
        _[:, sum(_, dims=1)[1, :] .!= 0] |> 
        _[sum(_, dims=2)[:, 1] .!= 0,:] 
    
    #Calculate beta diversity components
    components = beta_diversity(df_matrix; quant=quant)

    temporal_beta_div_summary = DataFrames.DataFrame(
        temporal_BDtotal =  components.BDtotal,
        temporal_Repl =  components.Repl,
        temporal_RichDif =  components.RichDif)

        return temporal_beta_div_summary
end

