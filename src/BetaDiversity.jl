#using Pkg
#Pkg.activate(".")

#using Distances
#using DataFrames
#using Statistics
#using Pipe: @pipe

###Beta Diveristy

"""
    beta_diversity(mat::Matrix; quant::Bool) -> DataFrame
Calculate beta diversity for given ecological data. This function supports both binary (presence/absence) and quantitative data. It returns a DataFrame containing the calculated beta diversity indices.

# Arguments
- `mat::Matrix`: A matrix where each row represents a sample and each column represents a species. The elements of the matrix should represent the presence/absence or abundance of species.
- `quant::Bool`: A boolean flag that indicates whether the data is quantitative. Default is `false`, which means the data is treated as binary.

# Returns
- `DataFrame`: A DataFrame with the following columns:
- `BDtotal`: Total beta diversity.
- `Repl`: Replacement component of diversity.
- `RichDif`: Richness difference component of diversity.

# Example
# Example
Load sample data included in the package and calculate beta diversity:
```julia
using CSV, DataFrames
using Pipe: @pipe
sample_data_path = joinpath(pkgdir(MetaCommunityMetrics), "data", "rodent_abundance_data.csv") # assign the path to the sample data
sample_matrix = @pipe CSV.read(sample_data_path, DataFrame; header=true) |> #read in the sample data
        filter(row -> row[:Sampling_date_order] == 1, _) |> # filter rows where Sampling_date_order == 1 and plot == 1
        select(_, Not(1:6)) |> # remove the first four columns
        Matrix(_) #covert the dataframe to a matrix
binary_result = beta_diversity(sample_matrix; quant=false) #calculate beta diversity for binary data
quant_result = beta_diversity(sample_matrix; quant=true) #calculate beta diversity for quantitative data

"""
#a function to calculate beta diversity, for binary/quantitative data
function beta_diversity(mat::Matrix; quant::Bool)
    #Error handling
    if isempty(mat)
        throw(ArgumentError("Input matrix is empty."))
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
#temporal mean of spatial beta-diversity,  beta-diversity among sites averaged across time

function mean_spatial_beta_div(abundance::Vector, time::Vector, patch::Vector, species::Vector; quant::Bool)
    #Create the required data frame
    df = DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=patch,
        Species=species)
    #Identify species that are present in all patches and time steps
    persist_sp_df=#idetify species with a total N >0 across  all patches all time steps. 
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

    for t in unique(dynamic_df.Time)
        #println("Time$t")
        subset_df=filter(row -> row[:Time]==t, dynamic_df)
        df_wide =@pipe unstack(subset_df, :Patch, :Species, :N) |> #pivot wider
        _[:,2:end] |> #Remove the patch column
        Matrix(_)#convert to matrix
        df_wide .= coalesce.(df_wide, 0) #replace missing values with zeros
        #Calculate beta diversity components
        components = beta_diversity(df_wide, quant)
        spatial_beta_div_df= [spatial_beta_div_df; components]
    end

    mean_spatial_beta_div_summary = DataFrames.DataFrame(
        mean_spatial_BDtotal = mean(spatial_beta_div_df.BDtotal),
        mean_spatial_Repl = mean(spatial_beta_div_df.Repl),
        mean_spatial_RichDif = mean(spatial_beta_div_df.RichDif))

    return mean_spatial_beta_div_summary
end
#spatial mean of temporal beta-diversity, beta-diversity along all time points averaged across site
function mean_temporal_beta_div(abundance::Vector, time::Vector, patch::Vector, species::Vector; quant::Bool)
    #Create the required data frame
    df = DataFrames.DataFrame(
        N=abundance,
        Time=time,
        Patch=patch,
        Species=species)
    #Identify species that are present in all patches and time steps
    persist_sp_df=#idetify species with a total N >0 across  all patches all time steps. 
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

    temporal_beta_div_df = DataFrames.DataFrame()  #a data frame to store beta diversity components at all patches

    for p in unique(dynamic_df.Patch)
        subset_df = filter(row -> row[:Patch]==p, dynamic_df)
        df_wide = @pipe unstack(subset_df, :Time, :Species, :N) |> #pivot wider
        _[:,2:end] |> #Remove the patch column
        Matrix(_)
        df_wide .= coalesce.(df_wide, 0)
         #Calculate beta diversity components
        components = beta_diversity(df_wide, quant)
        temporal_beta_div_df= [temporal_beta_div_df; components]
    end

    mean_temporal_beta_div_summary = DataFrames.DataFrame(
        mean_temporal_BDtotal = mean(temporal_beta_div_df.BDtotal),
        mean_temporal_Repl = mean(temporal_beta_div_df.Repl),
        mean_temporal_RichDif = mean(temporal_beta_div_df.RichDif))
    
    return mean_temporal_beta_div_summary
end

using Pkg
using Pipe: @pipe
Pkg.dev(path="/Users/jennycheung/.julia/dev/MetaCommunityMetrics")