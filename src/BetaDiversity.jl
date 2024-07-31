# src/BetaDiveristy.jl

"""
    beta_diversity(mat::Matrix; quant::Bool) -> DataFrame

Calculate beta diversity for a given biodiversity data. This function supports both binary (presence/absence) and quantitative data.
For binary data, the function calculates Podani family, Jaccard-based indices. For quantitative data, the function calculates Podani family, Ruzicka-based indices.
The function returns a DataFrame containing the calculated beta diversity indices. 
Empty patches have to be removed before calculation. Please refer to Example.jl for more details.

# Arguments
- `mat::Matrix`: A matrix where each row represents a sample and each column represents a species. The elements of the matrix should represent the presence/absence or abundance of species.
- `quant::Bool`: A boolean flag that indicates whether the data is quantitative. Default is `false`, which means the data is treated as binary.

# Returns
- `DataFrame`: A DataFrame with the following columns:
- `BDtotal`: Total beta diversity, which captures the overall dissimilarity between local communities.
- `Repl`: Replacement component of diversity, which reflects how many species are different in one site compared to another, ignoring the species that are mere additions or subtractions.
- `RichDif`: Richness difference component of diversity, which captures the disparity in biodiversity in terms of the count of species present, without taking into account the specific identities or distributions of those species.

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

Calculate the mean spatial beta diversity components of a metacommunity over time based on species abundances or presence-absences.

# Arguments
- `abundance::Vector`: A vector containing abundance data for each species across different samples.
- `time::Vector`: A vector indicating the time each sample was taken.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.
- `species::Vector`: A vector indicating the species associated with each abundance entry.
- `quant::Bool`: Optional boolean flag to indicate whether the data should be treated as quantitative (default is `false`, treating data as binary presence/absence).
    
# Returns
- `DataFrame`: A DataFrame containing the mean values of total beta diversity, replacement, and richness difference components across all time points. Columns are `mean_spatial_BDtotal`, `mean_spatial_Repl`, and `mean_spatial_RichDif`.
"""
function mean_spatial_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String};quant::Bool)
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

Calculate the mean temporal beta diversity components acorss all patches based on species abundances or presence-absences.

# Arguments
- `abundance::Vector`: A vector containing abundance data for each species across different samples.
- `time::Vector`: A vector indicating the time each sample was taken.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.
- `species::Vector`: A vector indicating the species associated with each abundance entry.
- `quant::Bool`: Optional boolean flag to indicate whether the data should be treated as quantitative (default is `false`, treating data as binary presence/absence).
    
# Returns
- `DataFrame`: A DataFrame containing the mean values of total beta diversity, replacement, and richness difference components across all pactches. Columns are `mean_temporal_BDtotal`, `mean_temporal_Repl`, and `mean_temporal_RichDif`.
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

