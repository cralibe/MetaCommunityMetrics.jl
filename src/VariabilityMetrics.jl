# src/VariabilityMetrics.jl 
###Community Dynamics metrics

using DataFrames, Pipe

"""
    CV_meta(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector) -> DataFrame

Calculates coefficients of variation (CV) for species and community biomass at both local and regional scales within a metacommunity.

Arguments
- `abundance::AbstractVector`: Vector representing the abundance of species.
- `time::AbstractVector`: Vector representing sampling dates.
- `site::AbstractVector`: Vector representing site names or IDs.
- `species::AbstractVector`: Vector representing species names or IDs.

Returns
- `DataFrame`: A DataFrame containing the following columns:
    - `CV_s_l`: Local-scale average species variability.
    - `CV_s_r`: Regional-scale average species variability.
    - `CV_c_l`: Local-scale average community variability.
    - `CV_c_r`: Regional-scale community variability.

Example
```@jildoctest
julia> using MetaCommunityMetrics, Pipe

julia> df = @pipe load_sample_data()
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

julia> CV_summary_df = CV_meta(df.Abundance, df.Sampling_date_order, df.plot, df.Species)
1×4 DataFrame
 Row │ CV_s_l   CV_s_r    CV_c_l    CV_c_r   
     │ Float64  Float64   Float64   Float64  
─────┼───────────────────────────────────────
   1 │ 1.48859  0.944937  0.718266  0.580183
```
"""
function CV_meta(abundance::AbstractVector, time::AbstractVector, site::AbstractVector, species::AbstractVector)

    ###Reorganize the data to only include the necessary columns
    df = DataFrames.DataFrame(
        Abundance=abundance,
        Time=time,
        Patch=site,
        Species=species)
    
    # Get unique values for each column
    species_ids = unique(df.Species)
    patches = unique(df.Patch)
    times = unique(df.Time)
    #Initialize matrices
    num_species = length(unique(df.Species))
    num_patches = length(unique(df.Patch))
    num_times = length(unique(df.Time))

    abundance_matrices = zeros(Float64, num_species, num_times, num_patches)

    # Fill the abundance matrices
    # Create mapping dictionaries for faster lookups
    species_to_idx = Dict(s => i for (i, s) in enumerate(species_ids))
    time_to_idx = Dict(t => i for (i, t) in enumerate(times))
    patch_to_idx = Dict(p => i for (i, p) in enumerate(patches))

    # Then use direct dictionary lookups
    for (idx, row) in enumerate(eachrow(df))
        i = species_to_idx[row.Species]
        j = time_to_idx[row.Time]
        k = patch_to_idx[row.Patch]
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
                sd_species_patch_ik[i,k] = std(abundance_matrices[i,:,k], corrected=true)
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