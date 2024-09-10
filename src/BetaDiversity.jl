# src/BetaDiveristy.jl

"""
    beta_diversity(mat::Matrix; quant::Bool) -> DataFrame

Calculate beta diversity for a given biodiversity data. This function supports both binary (presence/absence) and quantitative data.

Arguments
- `mat::Matrix`: A matrix where each row represents a sample and each column represents a species. The elements of the matrix should represent the presence/absence or abundance of species.
- `quant::Bool`: A boolean flag that specifies whether the data is quantitative. By default, it is set to `false`, which means the data will be treated as binary. In this case, any quantitative data will be converted to binary, and beta diversity is calculated using the Podani family’s Jaccard-based indices. If `true`, the data is treated as quantitative, and beta diversity is calculated using the Podani family’s Ruzicka-based indices. For binary data, `quant` must remain set to `false`.

Returns
- `DataFrame`: A DataFrame with the following columns:
    - `BDtotal`: Total beta diversity, which captures the overall dissimilarity between local communities.
    - `Repl`: Replacement component of diversity, which reflects how many species are different in one site compared to another, ignoring the species that are mere additions or subtractions.
    - `RichDif`: Richness difference component of diversity, which captures the disparity in biodiversity in terms of the count of species present, without taking into account the specific identities or distributions of those species.

Details
- Empty patches have to be removed before calculation.
- Species that were not recorded at the given time step have to be removed before calculation.
- For binary data, the function calculates Podani family, Jaccard-based indices. 
- For quantitative data, the function calculates Podani family, Ruzicka-based indices.




Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
48735×10 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64   
───────┼────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0
                                                                                          48725 rows omitted

julia> matrix_with_abundance =  @pipe df |>
           filter(row -> row[:Sampling_date_order] == 1, _) |> 
           select(_, Not(:Presence)) |>
           unstack(_, :Species, :Abundance, fill=0) |>
           select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> 
           Matrix(_) |> 
           x -> x[:, sum(x, dims=1)[1, :] .!= 0] |> 
           x -> x[vec(sum(x, dims=2)) .!= 0, :]
15×5 Matrix{Int64}:
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 2  0  0  0  0
 1  0  0  1  0
 4  0  0  1  0
 1  0  0  0  0
 0  1  0  0  1
 0  0  0  1  0
 0  0  0  2  0
 1  0  0  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0  1  0  0
 0  0  0  0  1

julia> matrix_with_presence = @pipe df |> 
           filter(row -> row[:Sampling_date_order] == 1, _) |> 
           select(_, Not(:Abundance)) |>
           unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
           select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> 
           Matrix(_) |>     
           x -> x[:, sum(x, dims=1)[1, :] .!= 0]
15×5 Matrix{Int64}:
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  1  0
 1  0  0  1  0
 1  0  0  0  0
 0  1  0  0  1
 0  0  0  1  0
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
    mean_spatial_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}; quant::Bool) -> DataFrame

Calculate the mean spatial beta diversity components of a metacommunity over time based on species abundances or presence-absences using the function `beta_diversity`.

Arguments
- `abundance::Vector`: A vector containing abundance data for each species across different samples.
- `time::Vector`: A vector indicating the time each sample was taken.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.
- `species::Vector`: A vector indicating the species associated with each abundance entry.
- `quant::Bool`: A boolean flag that specifies whether the data is quantitative. By default, it is set to `false`, which means the data will be treated as binary. In this case, any quantitative data will be converted to binary, and beta diversity is calculated using the Podani family’s Jaccard-based indices. If `true`, the data is treated as quantitative, and beta diversity is calculated using the Podani family’s Ruzicka-based indices. For binary data, `quant` must remain set to `false`.

Returns
- `DataFrame`: A DataFrame containing the mean values of total beta diversity, replacement, and richness difference components across all time points. Columns are `mean_spatial_BDtotal`, `mean_spatial_Repl`, and `mean_spatial_RichDif`.

Details
- This function uses the `beta_diversity` function to calculate beta diversity components for each time step.
- This function will remove the empty patches before calculating beta diversity.
- This function will remove species that were not recorded at the given time step before calculating beta diversity.
- For binary data, the function calculates Podani family, Jaccard-based indices. 
- For quantitative data, the function calculates Podani family, Ruzicka-based indices.

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
48735×10 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64   
───────┼────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0
                                                                                          48725 rows omitted

julia> result_using_abanduce_data_1 = mean_spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
1×3 DataFrame
 Row │ mean_spatial_BDtotal  mean_spatial_Repl  mean_spatial_RichDif 
     │ Float64               Float64            Float64              
─────┼───────────────────────────────────────────────────────────────
   1 │             0.353812           0.168584              0.185228
        
julia> result_using_abanduce_data_2 = mean_spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ mean_spatial_BDtotal  mean_spatial_Repl  mean_spatial_RichDif 
     │ Float64               Float64            Float64              
─────┼───────────────────────────────────────────────────────────────
   1 │             0.309073            0.16386              0.145214  

julia> result_using_binary_data = mean_spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ mean_spatial_BDtotal  mean_spatial_Repl  mean_spatial_RichDif 
     │ Float64               Float64            Float64              
─────┼───────────────────────────────────────────────────────────────
   1 │             0.309073            0.16386              0.145214
```

"""
function mean_spatial_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}; quant::Bool)
    #Create the required data frame
    df = DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=patch,
        Species=species)
    
    persist_sp_df=#idetify species with a total N >0 across all patches all time steps. 
    @pipe df |> #read in the model outputs from every landscape
    groupby(_, :Species) |> #group the dataframe by every species
    combine(_, :N => sum => :Total_N) |> #calculate the total abundance of every spesice within the whole sampling period
    filter(row -> row[:Total_N] > 0, _) #select rows that has a N>0.

    #Identify sites that are not empty at a time step
    non_empty_site = @pipe df |>
    groupby(_, [:Time, :Patch]) |>  # Group by time and patch
    combine(_, :N => sum => :Total_N) |>  # Sum N values within each time and patch combination
    filter(row -> row.Total_N > 0, _) 

    #Prepare a presence-absence matrix for beta diversity calculation
    dynamic_df = @pipe df |>
    innerjoin(persist_sp_df, _, on=:Species) |>  # Join with the original data to filter species
    select(_, Not(:Total_N)) |>  # Remove the total N column
    innerjoin(_, non_empty_site, on=[:Time, :Patch])|>  # Remove empty sites at one time step
    select(_, Not(:Total_N)) 


    spatial_beta_div_df = DataFrames.DataFrame()  #a data frame to store beta diversity components from evey time point
    for t in unique(df.Time)
        #println("Time$t")
        subset_df=filter(row -> row[:Time]==t, df)
        df_wide =@pipe unstack(subset_df, :Patch, :Species, :N) |> #pivot wider
        _[:,2:end] #Remove the patch column

         # Replace missing values with zeros before summing
         df_wide .= coalesce.(df_wide, 0)
 
         # Removing columns where the sum is zero
        non_zero_columns = names(df_wide)[.!iszero.(sum.(eachcol(df_wide)))]
        df_wide = select(df_wide, non_zero_columns)

         # Convert to matrix after filtering
        df_matrix = Matrix(df_wide)

        #Calculate beta diversity components
        components = beta_diversity(df_matrix; quant=quant)
        spatial_beta_div_df= [spatial_beta_div_df; components]
    end

    mean_spatial_beta_div_summary = DataFrames.DataFrame(
        mean_spatial_BDtotal = mean(spatial_beta_div_df.BDtotal),
        mean_spatial_Repl = mean(spatial_beta_div_df.Repl),
        mean_spatial_RichDif = mean(spatial_beta_div_df.RichDif))

    return mean_spatial_beta_div_summary
end

"""
    mean_temporal_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String};quant::Bool) -> DataFrame

Calculate the mean temporal beta diversity components acorss all patches based on species abundances or presence-absences using the function `beta_diversity`.

Arguments
- `abundance::Vector`: A vector containing abundance data for each species across different samples.
- `time::Vector`: A vector indicating the time each sample was taken.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.
- `species::Vector`: A vector indicating the species associated with each abundance entry.
- `quant::Bool`: A boolean flag that specifies whether the data is quantitative. By default, it is set to `false`, which means the data will be treated as binary. In this case, any quantitative data will be converted to binary, and beta diversity is calculated using the Podani family’s Jaccard-based indices. If `true`, the data is treated as quantitative, and beta diversity is calculated using the Podani family’s Ruzicka-based indices. For binary data, `quant` must remain set to `false`.

Returns
- `DataFrame`: A DataFrame containing the mean values of total beta diversity, replacement, and richness difference components across all pactches. Columns are `mean_temporal_BDtotal`, `mean_temporal_Repl`, and `mean_temporal_RichDif`.

Details
- This function uses the `beta_diversity` function to calculate beta diversity components for each patch.
- This function will remove the empty patches before calculating beta diversity.
- This function will remove species that were not recorded at the given time step before calculating beta diversity.
- For binary data, the function calculates Podani family, Jaccard-based indices. 
- For quantitative data, the function calculates Podani family, Ruzicka-based indices.

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames

julia> df = load_sample_data()
48735×10 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64   
───────┼────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0
                                                                                          48725 rows omitted

julia> result_using_abanduce_data_1 = mean_temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
1×3 DataFrame
 Row │ mean_temporal_BDtotal  mean_temporal_Repl  mean_temporal_RichDif 
     │ Float64                Float64             Float64               
─────┼──────────────────────────────────────────────────────────────────
   1 │               0.37186            0.152416               0.219444
        
julia> result_using_abanduce_data_2 = mean_temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ mean_temporal_BDtotal  mean_temporal_Repl  mean_temporal_RichDif 
     │ Float64                Float64             Float64               
─────┼──────────────────────────────────────────────────────────────────
   1 │              0.313244            0.139936               0.173309

julia> result_using_binary_data = mean_temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)
1×3 DataFrame
 Row │ mean_temporal_BDtotal  mean_temporal_Repl  mean_temporal_RichDif 
     │ Float64                Float64             Float64               
─────┼──────────────────────────────────────────────────────────────────
   1 │              0.313244            0.139936               0.173309
```

"""
function mean_temporal_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String};quant::Bool)
    #Create the required data frame
    df = DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=patch,
        Species=species)

    temporal_beta_div_df = DataFrames.DataFrame()  #a data frame to store beta diversity components at all patches

    for p in unique(df.Patch)
        subset_df = filter(row -> row[:Patch]==p, df)
        df_wide = @pipe unstack(subset_df, :Time, :Species, :N) |> #pivot wider
        _[:,2:end]  #Remove the patch column

         # Replace missing values with zeros before summing
        df_wide .= coalesce.(df_wide, 0)
 
         # Removing columns where the sum is zero
        non_zero_columns = names(df_wide)[.!iszero.(sum.(eachcol(df_wide)))]
        df_wide = select(df_wide, non_zero_columns)

         # Convert to matrix after filtering
        df_matrix = Matrix(df_wide)

         #Calculate beta diversity components
        components = beta_diversity(df_matrix, quant=quant)
        temporal_beta_div_df= [temporal_beta_div_df; components]
    end

    mean_temporal_beta_div_summary = DataFrames.DataFrame(
        mean_temporal_BDtotal = mean(temporal_beta_div_df.BDtotal),
        mean_temporal_Repl = mean(temporal_beta_div_df.Repl),
        mean_temporal_RichDif = mean(temporal_beta_div_df.RichDif))
    
    return mean_temporal_beta_div_summary
end

