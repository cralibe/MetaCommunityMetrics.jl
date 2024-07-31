using Pkg
Pkg.activate(".")

using Pipe: @pipe
using DataFrames

###Community Dynamics metrics

function CV_fun(abundance) # a function to calculated the general coffecieint of variation
    if (mean(abundance)) == 0
        #println("mean_abundance = 0")
        CV = 0
    else
        CV = std(abundance, corrected=true)/mean(abundance) #when corrected=true, the sum is scaled with n-1
    end
end

function CV_meta_compo_fun(dynamic_df)
    ##Variability metrics
    #Temporal variability of species i within the patch k
    species_patch_df = groupby(dynamic_df, [:Species,:Patch]) # Group by both Patch and Time
    species_patch_cv_df=combine(species_patch_df, :N => CV_fun => :CV)

    #Temporal variability of the metapopulation biomass of species i
    species_metacom_cv_df = 
    @pipe dynamic_df|> 
    groupby(_, [:Species, :Time]) |>
    combine(_, :N => sum => :meta_N) |>
    groupby(_, [:Species]) |>
    combine(_, :meta_N => CV_fun => :CV)

    #Temporal variability of total community biomass of patch k
    patch_cv_df = 
    @pipe dynamic_df |> 
    groupby(_, [:Patch, :Time]) |>
    combine(_, :N => sum => :patch_N) |>
    groupby(_, [:Patch]) |>
    combine(_, :patch_N => CV_fun => :CV)
    
    #Temporal mean biomass of the whole metacommunity
    meta_meanN_df = 
    @pipe dynamic_df |>
    groupby(_, [:Time]) |>
    combine(_, :N => sum => :meta_N)
    meta_meanN = mean(meta_meanN_df.meta_N)


    #Local-scale average species variability, defined as the weighted average of local population variability (CV_i_k) across species and patches
    species_patch_meanN_df = 
    @pipe species_patch_df |>
    combine(_, :N => mean => :mean_N) #Temporal mean biomass of species i in patch k
    species_patch_meanN_df.weight = species_patch_meanN_df.mean_N/meta_meanN #weights
    CV_s_l= mean(species_patch_cv_df.CV, weights(species_patch_meanN_df.weight))#average

    #Regional-scale average species variability,defined as the weighted average of metapopulation variability (CV_i_R) across species
    species_metacom_meanN_df =
    @pipe dynamic_df|> 
    groupby(_, [:Species, :Time]) |>
    combine(_, :N => sum => :meta_N) |>
    groupby(_, [:Species]) |> 
    combine(_, :meta_N => mean => :mean_N)#Temporal mean metapopulation biomass of species i
    species_metacom_meanN_df.weight=species_metacom_meanN_df.mean_N/meta_meanN #weights
    CV_s_r = mean(species_metacom_cv_df.CV, weights(species_metacom_meanN_df.weight)) #average

    #Local-scale average community variability,defined as the weighted average of community variability (CVC,k) across patches, the square of which corresponds to the alpha variability in Wang and Loreau (2014, 2016)
    patch_meanN_df =
    @pipe dynamic_df |> 
    groupby(_, [:Patch, :Time]) |>
    combine(_, :N => sum => :patch_N) |>
    groupby(_, [:Patch]) |>
    combine(_, :patch_N => mean => :mean_N) #Temporal mean community biomass of patch k
    patch_meanN_df.weight=patch_meanN_df.mean_N/meta_meanN
    CV_c_l = mean(patch_cv_df.CV, weights(patch_meanN_df.weight)) 


    #Regional-scale community variability or metacommunity variability, the square of which corresponds to the gamma variability in Wang and Loreau (2014, 2016)
    species_time_patch_df=
    @pipe dynamic_df |>
    select(_, :Time, :Patch, :Species, :N)|>
    sort(_, [:Time, :Species, :Patch])

    times = unique(species_time_patch_df.Time)
    species = unique(species_time_patch_df.Species)
    patch = unique(species_time_patch_df.Patch)

    T = length(times)
    S = length(species)
    P = length(patch)

    abundance_data = zeros(T, S, P)

    for row in eachrow(species_time_patch_df)
        t_index = findfirst(times .== row.Time)
        s_index = findfirst(species .== row.Species)
        p_index = findfirst(patch .== row.Patch)
        abundance_data[t_index, s_index, p_index] = row.N
    end

    # Initialize means
    mean_species_patch = mean(abundance_data[:, :, :], dims=1) #calculates the mean abundance for each species in each patch by averaging over time

    # Initialize covariance matrix
    Cov = zeros(Float64, S * P, S * P)

    # Calculate temporal covariance matrix using nested loops
    @showprogress 1 "Calculating temporal covariance matrix..." for i in 1:S
        for j in 1:S
            for k in 1:P
                for l in 1:P
                    for t in 1:T
                        Cov[(i - 1) * P + k, (j - 1) * P + l] +=
                        (abundance_data[t, i, k] - mean_species_patch[1, i, k]) * (abundance_data[t, j, l] - mean_species_patch[1, j, l])
                    end
                    # Normalize by (T-1)
                    if (i == j && k == l) || ((i - 1) * P + k > (j - 1) * P + l) #exclude the elements in the diagonal and below the diagonal of the covarience matrix.
                        Cov[(i - 1) * P + k, (j - 1) * P + l] = 0.0
                    else   
                        Cov[(i - 1) * P + k, (j - 1) * P + l] /= (T - 1)#assigning the division result back to the same element in the matrix.
                    end
                end
            end
        end
    end

    # Calculate the summation of Cov_(i,j,k,l) using nested loops
    #=total_sum = 0.0
    for i in 1:num_species
        for j in 1:num_species
            for k in 1:num_patches
                for l in 1:num_patches
                    total_sum += Cov[(i - 1) * num_patches + k, (j - 1) * num_patches + l]
                end
            end
        end
    end=#

    CV_c_r = sum(Cov)/meta_meanN

    CV_summary_df=DataFrames.DataFrame(
        CV_s_l=CV_s_l,
        CV_s_r=CV_s_r,
        CV_c_l=CV_c_l,
        CV_c_r=CV_c_r
    )
end

#####The simpler version (without covariences)
###Community Dynamics Metrics
function CV_meta_compo_fun(dynamic_df)

    # Get unique values for each column
    species_ids = unique(dynamic_df.Species)
    patches = unique(dynamic_df.Patch)
    times = unique(dynamic_df.Time)
    #Initialize matrices
    num_species = length(unique(dynamic_df.Species))
    num_patches = length(unique(dynamic_df.Patch))
    num_times = length(unique(dynamic_df.Time))

    abundance_matrices = zeros(Float64, num_species,num_times, num_patches)

    # Fill the abundance matrices
    for (idx, row) in enumerate(eachrow(dynamic_df))
        i = findfirst(species_ids .== row.Species)
        j = findfirst(times .== row.Time)
        k = findfirst(patches .== row.Patch)
        abundance_matrices[i, j, k] = row.N
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
end