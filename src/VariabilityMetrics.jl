# src/VariabilityMetrics.jl 
###Community Dynamics metrics

using DataFrames, Pipe

"""
    CV_meta(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}) -> DataFrame

Calculates various coefficients of variation (CV) for species and community biomass at both local and regional scales within a metacommunity.

Arguments
- `abundance::AbstractVector`: A vector representing the abundance of species.
- `time::AbstractVector`: A vector representing the time points at which the abundance measurements were taken.
- `patch::Union{AbstractVector, String}`: A vector or single value representing the patch or plot identifier.
- `species::Union{AbstractVector, String}`: A vector or single value representing the species identifier.

Returns
- `DataFrame`: A DataFrame containing the following columns:
    - `CV_s_l`: Local-scale average species variability.
    - `CV_s_r`: Regional-scale average species variability.
    - `CV_c_l`: Local-scale average community variability.
    - `CV_c_r`: Regional-scale community variability.

Details
This function calculates the coefficients of variation (CV) for species and community biomass at both local and regional scales. The calculation involves several steps:
1. **Reorganization of Data:** The input data is organized into a DataFrame with columns for abundance, time, plot, and species.
2. **Mean Calculations:** Temporal mean species abundance is calculated for each species in each patch, as well as the overall temporal mean biomass.
3. **Temporal Variance Calculations:** Temporal variance is calculated for each species within patches, for species across patches, for the community biomass within patches, and for the overall metacommunity biomass.
4. **CV Calculations:** The coefficients of variation are calculated for species and community biomass at both local and regional scales.
5. **Output:** The results are returned in a DataFrame summarizing the CVs for local and regional scales.

Example
```@jildoctest
julia> abundance = [10, 20, 15, 30, 25]
5-element Vector{Int64}:
 10
 20
 15
 30
 25

julia> time = [1, 1, 2, 2, 3]
5-element Vector{Int64}:
 1
 1
 2
 2
 3

julia> patch = ["A", "A", "A", "B", "B"]
5-element Vector{String}:
 "A"
 "A"
 "A"
 "B"
 "B"

julia> species = ["Sp1", "Sp2", "Sp1", "Sp2", "Sp1"]
5-element Vector{String}:
 "Sp1"
 "Sp2"
 "Sp1"
 "Sp2"
 "Sp1"

julia> CV_summary_df = CV_meta(abundance, time, patch, species)
1×4 DataFrame
 Row │ CV_s_l   CV_s_r    CV_c_l    CV_c_r   
     │ Float64  Float64   Float64   Float64  
─────┼───────────────────────────────────────
   1 │ 1.13137  0.870457  0.282843  0.223607

```
"""
function CV_meta(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String})
    ### Reorganize the data to only include the necessary columns
    df = DataFrame(
        Abundance = abundance,
        Time = time,
        plot = patch,
        Species = species
    )

    ### Mean calculations
    ## Temporal mean species abundance and number of unique time steps for each species in each patch
    temporal_mean_i_k_abundance = @pipe df |>
        unstack(_, :Species, :Abundance)|>
        coalesce.(_, 0) |>#Replace missing values with 0
        stack(_, Not([:plot, :Time]), variable_name = :Species, value_name = :Abundance) |>
        groupby(_, [:Species, :plot]) |>
        combine(_, :Abundance => mean => :mean_abundance)

    ## Temporal mean biomass of the whole
    temporal_mean_meta_meta_abundance = sum(temporal_mean_i_k_abundance.mean_abundance)

    ### Temporal variance calculations
    # A master dataframe to store the substraction results between the abundance and the mean abundance
    ab_minus_mean_ab_df = @pipe df |>
    unstack(_, :Species, :Abundance)|>
    coalesce.(_, 0) |>#Replace missing values with 0
    stack(_, Not([:plot, :Time]), variable_name = :Species, value_name = :Abundance) |>
    leftjoin(_, temporal_mean_i_k_abundance, on = [:Species, :plot]) |>
    transform(_, [:Abundance, :mean_abundance] => ((abundance, mean_abundance) -> (abundance .- mean_abundance)) => :ab_minus_mean_ab)

    # A master mutiplication dataframe
    master_multiplication_df = DataFrame(Time = Any[], Species_1 = Any[], Species_2 = Any[], plot_1 = Any[], plot_2 = Any[], multiplication = Float64[])

    for time in unique(ab_minus_mean_ab_df.Time)
        new_df = @pipe ab_minus_mean_ab_df |>
            filter(row -> row[:Time] == time, _) |>
            select(_, Not([:Time, :Abundance, :mean_abundance])) |>
            unstack(_, :Species, :ab_minus_mean_ab)

        # Get species columns (excluding the first column which is the plot identifier)
        sp_columns = names(new_df)[2:end]

        # Create combinations of species columns
        sp_combinations = [(sp1, sp2) for sp1 in sp_columns for sp2 in sp_columns]

        # Loop through each combination of species
        for (sp1, sp2) in sp_combinations
            sp1_rows = new_df[:, sp1]
            sp2_rows = new_df[:, sp2]

            # Generate all possible pairs of values from sp1_rows and sp2_rows and multiply them
            for i in 1:length(sp1_rows)
                for j in 1:length(sp2_rows)
                    plot_1 = new_df.plot[i]
                    plot_2 = new_df.plot[j]

                    plot1_value = sp1_rows[i]
                    plot2_value = sp2_rows[j]

                    # Calculate v_meta_meta for the current combination
                    multiplication = plot1_value * plot2_value

                    # Create a row for the current combination
                    row = DataFrame(Time = [time], Species_1 = [sp1], Species_2 = [sp2], plot_1 = [plot_1], plot_2 = [plot_2], multiplication = [multiplication])

                    # Append the row to the final DataFrame
                    append!(master_multiplication_df, row, promote = true)
                end
            end
        end
    end

    ## Temporal variance of species i in patch k
    v_ii_kk_df = @pipe master_multiplication_df |>
        filter(row -> row.Species_1 == row.Species_2, _)|>
        filter(row -> row.plot_1 == row.plot_2, _)|>
        groupby(_, [:Species_1, :Species_2, :plot_1, :plot_2]) |>
        combine(_, [:multiplication, :Time] => ((multiplication, Time) -> sum(multiplication) / (length(unique(Time)) - 1)) => :v_ii_kk) |>
        transform(_, :v_ii_kk => (ByRow(x -> (isnan(x) || isinf(x)) ? missing : x)) => :v_ii_kk) |> #NaN/Inf appers when there is only one time step, replacing them with missing
        dropmissing(_) #missing value is removed

    ## Temporal variance of metapopulation biomass (all plots) for species i
    v_ii_meta_df = @pipe master_multiplication_df |>
        filter(row -> row.Species_1 == row.Species_2, _) |>
        groupby(_, [:Species_1, :Species_2, :plot_1, :plot_2]) |>
        combine(_, [:multiplication, :Time] => ((multiplication, Time) -> sum(multiplication) / (length(unique(Time)) - 1)) => :v_ii_kl) |>
        transform(_, :v_ii_kl => (ByRow(x -> (isnan(x) || isinf(x)) ? missing : x)) => :v_ii_kl) |> #NaN/Inf appers when there is only one time step, replacing them with missing
        dropmissing(_) |> #missing value is removed
        groupby(_, [:Species_1, :Species_2]) |>
        combine(_, [:v_ii_kl] => sum => :v_ii_meta)
 

    ## Temporal variance of total community biomass of patches k 
    v_meta_kk_df = @pipe master_multiplication_df |>
        filter(row -> row.plot_1 == row.plot_2, _)|>
        groupby(_, [:Species_1, :Species_2, :plot_1, :plot_2]) |>
        combine(_, [:multiplication, :Time] => ((multiplication, Time) -> sum(multiplication) / (length(unique(Time)) - 1)) => :v_ij_kk) |>
        transform(_, :v_ij_kk => (ByRow(x -> (isnan(x) || isinf(x)) ? missing : x)) => :v_ij_kk) |> #NaN/Inf appers when there is only one time step, replacing them with missing
        dropmissing(_) |> #missing value is removed
        groupby(_, [:plot_1]) |>
        combine(_, [:v_ij_kk] => sum => :v_meta_kk)

    ## Temporal variance of the whole metacommunity
    v_meta_meta_df = @pipe master_multiplication_df |>
    groupby(_, [:Species_1, :Species_2, :plot_1, :plot_2]) |>
    combine(_, [:multiplication, :Time] => ((multiplication, Time) -> sum(multiplication) / (length(unique(Time)) - 1)) => :v_ij_kl) |>
    transform(_, :v_ij_kl => (ByRow(x -> (isnan(x) || isinf(x)) ? missing : x)) => :v_ij_kl) |>
    dropmissing(_) |>
    combine(_, [:v_ij_kl] => sum => :v_meta_meta)

    ### Local-scale average species variability
    CV_s_l = sum(sqrt.(v_ii_kk_df.v_ii_kk)) / temporal_mean_meta_meta_abundance
    ### Regional-scale average species variability
    CV_s_r = sum(sqrt.(v_ii_meta_df.v_ii_meta)) / temporal_mean_meta_meta_abundance
    ### Local-scale average community variability
    CV_c_l = sum(sqrt.(v_meta_kk_df.v_meta_kk)) / temporal_mean_meta_meta_abundance
    ### Regional-scale community variability
    CV_c_r = sqrt(sum(v_meta_meta_df.v_meta_meta)) / temporal_mean_meta_meta_abundance

    CV_summary_df = DataFrame(
        CV_s_l = CV_s_l,
        CV_s_r = CV_s_r,
        CV_c_l = CV_c_l,
        CV_c_r = CV_c_r
    )
    return CV_summary_df
end

"""
    CV_meta_simple(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}) -> DataFrame

Calculates coefficients of variation (CV) for species and community biomass at both local and regional scales within a metacommunity, using a simpler approach optimized for handling larger datasets.

Arguments
- `abundance::AbstractVector`: A vector representing the abundance of species.
- `time::AbstractVector`: A vector representing the time points at which the abundance measurements were taken.
- `patch::Union{AbstractVector, String}`: A vector or single value representing the patch or plot identifier.
- `species::Union{AbstractVector, String}`: A vector or single value representing the species identifier.

Returns
- `DataFrame`: A DataFrame containing the following columns:
    - `CV_s_l`: Local-scale average species variability.
    - `CV_s_r`: Regional-scale average species variability.
    - `CV_c_l`: Local-scale average community variability.
    - `CV_c_r`: Regional-scale community variability.

Details
This function is a simplified version of the `CV_meta` function, designed to efficiently handle larger datasets by avoiding complex covariance calculations. The steps include:

1. **Reorganization of Data:** The input data is organized into a DataFrame with columns for abundance, time, plot, and species, and then transformed into a 3D abundance matrix.
2. **Total Abundance Calculations:** The function calculates total abundances for species across time, within each patch, and for the entire metacommunity.
3. **Standard Deviation (SD) Calculations:** Temporal standard deviations of abundance are computed for the entire metacommunity, each patch, and each species.
4. **CV Calculations:** The coefficients of variation are calculated for species and community biomass at both local and regional scales.
5. **Output:** The results are returned in a DataFrame summarizing the CVs for local and regional scales.

Example
```@jildoctest
julia> abundance = [10, 20, 15, 30, 25]
5-element Vector{Int64}:
 10
 20
 15
 30
 25

julia> time = [1, 1, 2, 2, 3]
5-element Vector{Int64}:
 1
 1
 2
 2
 3

julia> patch = ["A", "A", "A", "B", "B"]
5-element Vector{String}:
 "A"
 "A"
 "A"
 "B"
 "B"

julia> species = ["Sp1", "Sp2", "Sp1", "Sp2", "Sp1"]
5-element Vector{String}:
 "Sp1"
 "Sp2"
 "Sp1"
 "Sp2"
 "Sp1"

julia> CV_summary_df = CV_meta_simple(abundance, time, patch, species)
1×4 DataFrame
 Row │ CV_s_l   CV_s_r    CV_c_l    CV_c_r  
     │ Float64  Float64   Float64   Float64 
─────┼──────────────────────────────────────
   1 │ 1.52817  0.687386  0.932183  0.31225
   
```
"""
function CV_meta_simple(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String})

    ###Reorganize the data to only include the necessary columns
    metacomm_df = DataFrames.DataFrame(
        Abundance=abundance,
        Time=time,
        Patch=patch,
        Species=species)
    
    # Get unique values for each column
    species_ids = unique(metacomm_df.Species)
    patches = unique(metacomm_df.Patch)
    times = unique(metacomm_df.Time)
    #Initialize matrices
    num_species = length(unique(metacomm_df.Species))
    num_patches = length(unique(metacomm_df.Patch))
    num_times = length(unique(metacomm_df.Time))

    abundance_matrices = zeros(Float64, num_species,num_times, num_patches)

    # Fill the abundance matrices
    for (idx, row) in enumerate(eachrow(metacomm_df))
        i = findfirst(species_ids .== row.Species)
        j = findfirst(times .== row.Time)
        k = findfirst(patches .== row.Patch)
        abundance_matrices[i, j, k] = row.Abundance
    end

    #total abundance of all species in the same time point
        # Initialize a vector to store total abundance of all species in the same time point
        ts_metacom = zeros(num_times)
        # Calculate the total abundance
        for t in 1:num_times
            ts_metacom[t] = sum(abundance_matrices[:, t, :])
        end

    # total abundance of all species in the same patch at the same time point
        # Initialize a vector to store the total abundance of all species in the same patch at the same time point
        ts_patch = zeros(Float64,num_times,num_patches)
        # Calculate the total abundance
        for t in 1:num_times
            for k in 1:num_patches
                ts_patch[t,k] = sum(abundance_matrices[:, t, k])
            end
        end
        
    # total abundance of species i at the same time point 
        # Initialize a vector to store the total abundance of total abundance of species i at the same time point 
        ts_species = zeros(Float64,num_species,num_times)
        # Calculate the total abundance
        for i in 1:num_species
            for t in 1:num_times
                ts_species[i,t] = sum(abundance_matrices[i, t, :])
            end
        end
    
    sd_metacom = std(ts_metacom, corrected=true)#sd of metacommunity across time
    sd_patch_k = [std(ts_patch[:, k], corrected=true) for k in 1:num_patches]#sd of total abundance across time for every patch
    sd_species_i = [std(ts_species[i, :], corrected=true) for i in 1:num_species]#sd of total number of species i across time
    
    # sd of species i across time for every patch
        #Initialize a vector to store the sd of total abundance of species i at the same time point  
        sd_species_patch_ik = zeros(Float64,num_species,num_patches)
        #Calculate the sd
        for i in 1:num_species
            for k in 1:num_patches
                sd_species_patch_ik[i,k] = std(abundance_matrices[i,:,k])
            end
        end 
      
    #metacommunity abundance averaged across time
    mean_metacom = mean(ts_metacom)
    
    CV_s_l = sum(sd_species_patch_ik)/mean_metacom
    CV_c_l = sum(sd_patch_k)/mean_metacom
    CV_s_r = sum(sd_species_i)/mean_metacom
    CV_c_r = sd_metacom/mean_metacom
  

    CV_summary_df=DataFrames.DataFrame(
        CV_s_l=CV_s_l,
        CV_s_r=CV_s_r,
        CV_c_l=CV_c_l,
        CV_c_r=CV_c_r
    )
    return CV_summary_df
end