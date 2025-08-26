# src/Internal.jl

module Internal

using DataFrames
using Random
using Distances
using LinearAlgebra
using Combinatorics
using Statistics
using StatsBase
using Pipe

### Internal function
## These functions is for internal use and supports the public functions.

## Internal fuctions for DNCI.jl
# A function assign single site to the nearest group
function assign_single_site_to_nearest_group(grouped_df, single_site, min_group)
    # Extract latitude and longitude of the single site
    single_lat = single_site.Latitude[1]
    single_lon = single_site.Longitude[1]
    
    min_distance = Inf
    nearest_group_id = nothing

    # Filter for other groups
    other_group = grouped_df[grouped_df.Group .!= min_group, :]

    
    other_group_with_xy = @pipe other_group |>
                            select(_, [:Site, :Latitude, :Longitude, :Group]) |>
                            unique(_)

    # Iterate over other groups
    for row in eachrow(other_group_with_xy)
        dist = sqrt((single_lat - row.Latitude)^2 + (single_lon - row.Longitude)^2)
        if dist < min_distance
            min_distance = dist
            nearest_group_id = row.Group
        end
    end

    if isnothing(nearest_group_id)
        error("No valid nearest group found for the single site.")
    end

    # Assign the single site to the nearest group
    #single_site.Group[1] = nearest_group_id

    return nearest_group_id
end

# A function that find the nearest site to the group with the fewest sites
function find_nearest_site(grouped_df, min_group)
    # Extract coordinates of sites in the group with the fewest species
    sites_in_min_group = grouped_df[grouped_df.Group .== min_group, :]
    
    # Check if there are any sites in the minimum group
    if isempty(sites_in_min_group)
        error("No sites found in the minimum group.")
    end

    # Calculate centroid of the cluster
    centroid_latitude = mean(unique(sites_in_min_group.Latitude))
    centroid_longitude = mean(unique(sites_in_min_group.Longitude))
    
    # Initialize variables to keep track of the nearest site and distance
    min_distance = Inf
    site_with_minimum_distance = nothing  # Initialize with nothing to handle cases where no site is found

    # Filter for other groups
    other_group = grouped_df[grouped_df.Group .!= min_group, :]

    # Check if there are other groups available
    if isempty(other_group)
        error("No other groups are avaliable to compare.")
    end

    other_group_with_xy = @pipe other_group |>
    select(_, [:Site, :Latitude, :Longitude, :Group]) |>
    unique(_)


    # Iterate over each site in other_group
    for row in eachrow(other_group_with_xy)
        # Calculate the Euclidean distance between the centroid and the current site
        dist = sqrt((centroid_latitude - row.Latitude)^2 + (centroid_longitude - row.Longitude)^2)
     
        # Update the nearest site if this site is closer
        if dist < min_distance
            min_distance = dist
            site_with_minimum_distance = row
        end
    end

    # Check if a nearest site was found
    if isnothing(site_with_minimum_distance)
        error("No nearest site could be determined.")
    end
        
    return site_with_minimum_distance
end
# A function that checks the condition for the groupings only
function check_conditions(grouped_df::DataFrame)
    sites_per_group_df = @pipe grouped_df |>
    groupby(_, :Group) |>
    combine(_, :Site => (x -> length(unique(x))) => :sites_per_group)

    # Calculate total richness per group
    total_richness_per_group = @pipe grouped_df |>
                groupby(_, [:Species, :Group]) |>
                combine(_,:Presence=>sum=>:Total_Presence)|>
                transform(_, :Total_Presence => ByRow(x->ifelse.(x> 0, 1, 0)) => :Presence)|> 
                groupby(_, :Group) |>
                combine(_, :Presence => sum => :Total_Richness)

    unique_groups = unique(grouped_df.Group)

    # Create lookup dictionaries for consistent ordering
    sites_dict = Dict(zip(sites_per_group_df.Group, sites_per_group_df.sites_per_group))
    richness_dict = Dict(zip(total_richness_per_group.Group, total_richness_per_group.Total_Richness))

    # Extract values as vectors in consistent order
    sites_vals = [sites_dict[g] for g in unique_groups]
    richness_vals = [richness_dict[g] for g in unique_groups]

    # Create pairwise difference matrices
    sites_diffs = abs.(sites_vals .- sites_vals')
    richness_diffs = abs.(richness_vals .- richness_vals')

    # Create max matrices for normalization
    sites_maxes = max.(sites_vals, sites_vals')
    richness_maxes = max.(richness_vals, richness_vals')

    # Calculate ratio matrices
    sites_ratios = sites_diffs ./ sites_maxes
    richness_ratios = richness_diffs ./ richness_maxes

    # Check condition (exclude diagonal since comparing group with itself gives 0)
    mask = .!(I(length(unique_groups)))  # Boolean mask excluding diagonal
    return (all((sites_ratios[mask] .<= 0.3) .& (richness_ratios[mask] .<= 0.4)))

    if any(value < 5 for value in values(sites_dict)) || length(unique_groups) <= 1
        return false
    end

    return true
end

# A function that checks the condition for the groupings and fix the groupings
function check_condition_and_fix(grouped_df)
    max_iterations = 1000  # Define a maximum number of iterations to prevent infinite loops
    iteration = 0

    while iteration < max_iterations

        # Group by 'Group' and calculate the number of unique sites per group
        sites_per_group_df = @pipe grouped_df |>
        groupby(_, :Group) |>
        combine(_, :Site => (x -> length(unique(x))) => :sites_per_group)

        # Calculate total richness per group
        total_richness_per_group = @pipe grouped_df |>
                    groupby(_, [:Species, :Group]) |>
                    combine(_,:Presence=>sum=>:Total_Presence)|>
                    transform(_, :Total_Presence => ByRow(x->ifelse.(x> 0, 1, 0)) => :Presence)|> 
                    groupby(_, :Group) |>
                    combine(_, :Presence => sum => :Total_Richness)

        unique_groups = unique(grouped_df.Group)

        # Create lookup dictionaries for consistent ordering
        sites_dict = Dict(zip(sites_per_group_df.Group, sites_per_group_df.sites_per_group))
        richness_dict = Dict(zip(total_richness_per_group.Group, total_richness_per_group.Total_Richness))

        # Extract values as vectors in consistent order
        sites_vals = [sites_dict[g] for g in unique_groups]
        richness_vals = [richness_dict[g] for g in unique_groups]

        # Create pairwise difference matrices
        sites_diffs = abs.(sites_vals .- sites_vals')
        richness_diffs = abs.(richness_vals .- richness_vals')

        # Create max matrices for normalization
        sites_maxes = max.(sites_vals, sites_vals')
        richness_maxes = max.(richness_vals, richness_vals')

        # Calculate ratio matrices
        sites_ratios = sites_diffs ./ sites_maxes
        richness_ratios = richness_diffs ./ richness_maxes

        # Check condition (exclude diagonal since comparing group with itself gives 0)
        mask = .!(I(length(unique_groups)))  # Boolean mask excluding diagonal
        condition_met = all((sites_ratios[mask] .<= 0.3) .& (richness_ratios[mask] .<= 0.4))


        if any(value < 5 for value in values(sites_dict)) || length(unique_groups) <= 1
            condition_met = false
        end
        
        if condition_met
            #println("Condition met, stopping iteration.")
            break
        else
             
            #println("Condition not met, adjusting data...") 

            min_group = first(argmin(sites_dict))

            # Check if there is only one site in the minimum group
            sites_in_min_group = grouped_df[grouped_df.Group .== min_group, :]
            if length(unique(sites_in_min_group.Site)) == 1
            # Find the nearest group to this single site
            grouped_df.Group[grouped_df.Group .== min_group, :] .= assign_single_site_to_nearest_group(grouped_df, sites_in_min_group, min_group)
            else
            site_with_minimum_distance = find_nearest_site(grouped_df, min_group) #Find the nearest site to the group with the fewest sites
            grouped_df.Group[grouped_df.Site .== site_with_minimum_distance.Site] .= min_group #assign the nearest site to the group with the fewest sites
            end
        end

        iteration += 1
    end

    #=if iteration == max_iterations
        println("Reached maximum iterations without meeting condition.")
    end=#

    return grouped_df
end

# A function to generate random binary matrices acoording to the assigned algorithm
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function permatfull(data::Matrix{Int}, fixedmar::String, Nperm::Int) #only for present/absence data

    #Assign for different permutation algorithms
    ALGO = ""
    if fixedmar == "rows"
        ALGO = "r0" #non-sequential algorithm for binary matrices 
        #that preserves the site (row) frequencies.
    elseif fixedmar == "columns"
        ALGO= "c0" #non-sequential algorithm for binary matrices 
        #that preserves species frequencies (Jonsson 2001).
    elseif fixedmar == "both"
        ALGO = "quasiswap" #non-sequential algorithm for binary matrices 
        #that implements a method where matrix is first filled honouring 
        #row and column totals, but with integers that may be larger than one. 
        #Then the method inspects random 2*2 matrices and performs a quasiswap on them.
    else
        throw(ArgumentError("Invalid fixedmar value, can only be `rows`, `columns` or `both`"))
    end
    
    permat=nullmodel(data, ALGO, Nperm)
    
    if !check_sums(data,permat, ALGO)
        if ALGO == "r0"
            throw(ArgumentError("The row sums do not match with the original matrix."))
        elseif ALGO == "c0"
            throw(ArgumentError("The column sums do not match with the original matrix."))
        elseif ALGO == "both"
            throw(ArgumentError("The row and column sums do not match with the original matrix."))
        end
    end
    
    return permat
    
end


    
# A function to generate null models acoording to the assigned algorithm
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function nullmodel(m::AbstractMatrix{Int}, ALGO::AbstractString, times::Int)
    if ALGO == "r0"
        return nullmodel_r0(m, times)
    elseif ALGO == "c0"
        return nullmodel_c0(m, times)
    elseif ALGO == "quasiswap"
        return nullmodel_quasiswap(m, times)
    else
        throw(ArgumentError("Invalid ALGO value"))
    end
end

# A function to generate random binary matrices with preserved row sums
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function nullmodel_r0(comm::AbstractMatrix{Int}, times::Int)
    
    n, m = size(comm)
    nullmodels = Vector{Matrix{Int}}(undef, times)

    for t in 1:times
    # Calculate row sums
    row_sums = sum(comm, dims=2)
    # Generate random binary matrix with preserved row sums
    nullmodel = zeros(Int, n, m)
        for i in 1:n
            nullmodel[i, shuffle(1:m)[1:row_sums[i]]] .= 1
        end
            nullmodels[t] = nullmodel
        end
    
    return nullmodels
end
# A function to generate random binary matrices with preserved column sums
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function nullmodel_c0(comm::AbstractMatrix{Int}, times::Int)

    n, m = size(comm)
    nullmodels = Vector{Matrix{Int}}(undef, times)
    
    for t in 1:times
        # Calculate column sums
        col_sums = sum(comm, dims=1)
        # Generate random binary matrix with preserved column sums
        nullmodel = zeros(Int, n, m)
        for j in 1:m
            nullmodel[shuffle(1:n)[1:col_sums[j]], j] .= 1
        end
        nullmodels[t] = nullmodel
    end
    
    return nullmodels
end
# A function to generate random binary matrices with preserved row and column sums using the Quasiswap algorithm
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function nullmodel_quasiswap(comm::AbstractMatrix{Int}, times::Int)

    n, m = size(comm)
    nullmodels = Vector{Matrix{Int}}(undef, times)

    for t in 1:times
        # Initialize nullmodel with the original matrix to preserve sums.
        nullmodel = copy(comm)
        
        # Randomly select pairs for potential swaps multiple times.
        for _ in 1:(10 * n * m) # Increase iterations for more thorough randomization.
            # Select two distinct rows and columns at random.
            rows = randperm(n)[1:2]
            cols = randperm(m)[1:2]
            
            # Extract the 2x2 submatrix.
            submatrix = nullmodel[rows, cols]
            
            # Determine if a quasiswap is possible while maintaining row/col sums.
            if sum(submatrix) == 2 && (submatrix[1, 1] + submatrix[2, 2] == 2 || submatrix[1, 2] + submatrix[2, 1] == 2)
                # Perform quasiswap if it leads to a valid configuration.
                nullmodel[rows, cols] .= 1 .- submatrix
            end
        end

        nullmodels[t] = nullmodel
    end
    
    return nullmodels
end

# A function to compare the original matrix with the null models in terms of row and column sums
function check_sums(original::AbstractMatrix{Int}, nullmodels::Vector{Matrix{Int}}, ALGO::AbstractString)
    original_row_sums = sum(original, dims=2)
    original_col_sums = sum(original, dims=1)

    for nullmodel in nullmodels
        if ALGO == "r0"
            nullmodel_row_sums = sum(nullmodel, dims=2)
            if !all(nullmodel_row_sums .== original_row_sums)
                return false
            end
        elseif ALGO == "c0"
            nullmodel_col_sums = sum(nullmodel, dims=1)
            if !all(nullmodel_col_sums .== original_col_sums)
                return false
            end
        elseif ALGO == "quasiswap"
            nullmodel_row_sums = sum(nullmodel, dims=2)
            nullmodel_col_sums = sum(nullmodel, dims=1)
            if !all(nullmodel_row_sums .== original_row_sums) || !all(nullmodel_col_sums .== original_col_sums)
                return false
            end
        else
            throw(ArgumentError("Invalid ALGO value: $ALGO"))
        end
    end

    return true
end

# A function to calculate Species contribution to average between-group dissimilarity.
# This function is a translation/adaptation of a function from the R package `vegan`, licensed under GPL-2 or later.
# Original package and documentation available at: https://cran.r-project.org/web/packages/vegan/index.html
function simper(comm::Matrix, groups::Vector)
    # Set EPS to square root of machine epsilon
    EPS = sqrt(eps(Float64))

    # Here is different from the orginal simper() in R
    # We use Zero-adjusted Bray-Curtis instead of the original Bray-Curtis
    # We added pseudo-species with value 1 to all sites 
    val = 1
    
    # Add pseudo-species column to every site
    comm_with_pseudo_species = hcat(comm, fill(val, size(comm, 1)))

    # Create a lower triangular matrix indicating whether each element (i, j) satisfies i > j
    tri = [i > j for i in 1:size(comm_with_pseudo_species, 1), j in 1:size(comm_with_pseudo_species, 1)]
    
    ## Species contributions of differences needed for every species,
    ## but denominator is constant. Bray-Curtis is actually
    ## manhattan/(mean(rowsums)) and this is the way we collect data
    # Calculate row sums to obtain the total diversity of the whole community at each site
    rs = sum(comm_with_pseudo_species, dims=2)
    # Calculate pairwise sums and extract lower triangular part
    pairwise_sums = rs .+ transpose(rs)
    result = pairwise_sums[tri]

    # Initialize an empty array to store species contributions
    spcontr = Matrix{Float64}(undef, 0, 0)
    
    # Iterate over each pair of sites in the community matrix
    for col_index in 1:size(comm_with_pseudo_species, 2)
        # Extract the ith column of the community matrix
        column_i = comm_with_pseudo_species[:, col_index:col_index]
        ZAP = 1e-15 #a threshold for setting small distances to zero in the distance matrix
        n_rows = size(column_i, 1)
        d = zeros(n_rows, n_rows)  # Initialize distance matrix
         # Calculate pairwise distances only for the lower triangular part
        for i in 2:n_rows
            for j in 1:i-1
                # Compute Manhattan distance between rows i and j
                distance = sum(abs.(column_i[i] .- column_i[j]))
                d[i, j] = distance < ZAP ? 0 : distance
            end
        end
        # Extract the lower triangle of the distance matrix and store it as a vector
        lower_triangle = d[tril(trues(size(d)), -1)]
        lower_triangle = reshape(lower_triangle, length(lower_triangle), 1)
        # Concatenate the lower triangle vector to spcontr
        if isempty(spcontr)
            spcontr = lower_triangle
        else
            spcontr = hcat(spcontr, lower_triangle)
        end
    end
    # Divide every value in each column of spcontr by the corresponding element in result
    spcontr ./= result

    #remove the last column of spcontr, which is the pseudo-species
    spcontr = spcontr[:, 1:end-1]

    #Get all combinations of 2 elements from unique_group
    comp = collect(combinations(unique(groups), 2))

    outlist = Dict{String, Any}() # Initialize an empty dictionary to store the output


    # function to match constrasts
    function contrmatch(X, Y, patt)
                
        return (X != Y) && any(in(X, p) for p in patt) && any(in(Y, p) for p in patt)

    end
   
    # Initialize an empty vector to store averages
    average_values = Matrix{Float64}[]
    # Iterate over patterns
    for k in 1:size(comp, 1)
        # Extract the current pattern
        patt = comp[k, :]
        # Initialize tmat as an array of booleans
        tmat = falses(length(groups), length(groups))
        # Compute tmat for the current pattern
        for i in 1:length(groups)
            for j in 1:length(groups)
                tmat[i, j] = contrmatch(groups[i], groups[j], patt)
            end
        end
        take = tmat[tril(trues(size(tmat)), -1)]
        # Initialize an array to store the selected rows for each column
        selected_rows = Vector[]

        # Iterate over each column of spcontr
        for col in eachcol(spcontr)
            # Extract the selected rows based on the true values in take
            selected = col[take]
            push!(selected_rows, selected)
        end
        selected_matrix = hcat(selected_rows...)
        
        #Species contribution to average between-group dissimilarity
        average=mean(selected_matrix, dims=1)

        # Store the average for the current k
        push!(average_values, average)
    end
    return average_values
end

# A function to correct dp4 matrix in the PerSIMPER function
function correct_permutation(comm::Matrix)
    SWAPcount = 1
    v = true
    dp4 = permatfull(comm, "columns",1) 

    while any(sum(dp4[1], dims=2) .== 0)#if any row sum = 0
        SWAPcount += 1
        v = false
        dp4 = permatfull(comm, "columns",1)

        if v == false && SWAPcount > 200
            while true
                for j in 1:size(dp4[1], 1)
                    if sum(dp4[1][j, :]) == 0 #if any row sum= 0
                        v = false
                        # Find species with the largest dispersal capacity (last quartile)
                        colSums = sum(dp4[1], dims=1)
                        highDispersalCols = findall(vec(colSums) .>= quantile(vec(colSums), 0.75)) #this indicates the position of the species with the largest dispersal capacity
                        tempHighCol = sample(highDispersalCols, 1, replace=false)[1] #a single random element from the array highDispersalCols without the possibility of choosing the same element again
                        # Find rich localities where the species with the largest dispersal capacity is present
                        presentCells = findall(dp4[1][:, tempHighCol] .> 0) #this indicates the position of the sites where the species with the largest dispersal capacity is present
                        richLocalities = 
                        presentCells[findall(sum(dp4[1][presentCells, :], dims=2) .>= 
                        quantile(vec(sum(dp4[1][presentCells, :], dims=2)), 0.75))]

                        selectedCell = sample(richLocalities, 1, replace=false)[1]

                        # Perform the swap
                        dp4[1][selectedCell, tempHighCol] = 0
                        dp4[1][j, tempHighCol] = 1

                        if all(sum(dp4[1], dims=1) .!= 0)
                            v = true
                            break
                        end
                    end
                end
                if v
                    break
                end
            end
        end
        
        if all(sum(dp4[1], dims=1) .!= 0)
            break
        end
    end

    return dp4
end

#PerSIMPER function : identification of the main assembly process; adpated from Corentin Gibert
#This function is a translation/adaptation of a function from the R package `DNCImper`, licensed under GPL-3.
#Original package and documentation available at: https://github.com/Corentin-Gibert-Paleontology/DNCImper
function PerSIMPER(comm::Matrix, groups::Vector, Nperm::Int=1000; count::Bool=false) #only for present/absence data

    AnaSimp = simper(comm, groups)#Species contribution to average between-group dissimilarity on the compared groups

    Contribution = sort(vec(AnaSimp[1]), rev=true)#Replication in a vector (named 'Contribution') of the sorting 
    #of species by their contribution to overall dissimilarity (OAD)

    Pourcent_Contribution = ((Contribution)/sum(Contribution))*100 #Conversion as a percentage of each species' contribution to the OAD

    #Randomization of the original community matrix 
    dp2 = permatfull(comm, "both", Nperm)
    dp3 = permatfull(comm, "rows", Nperm)
    dp4 = permatfull(comm, "columns", Nperm)
    
    #Generating matrices that will store the results (the ranked contribution of species to the OAD)
    #of the 1000 permutations of the original community matrix
    df2 = zeros(Nperm, size(comm, 2))
    df3 = zeros(Nperm, size(comm, 2))
    df4 = zeros(Nperm, size(comm, 2))

    for i in 1:Nperm
        
        if count == true && (i < 100 || i > round(Int, Nperm * 0.9))
            println(i)
        end
        
        #local dp4_filtered, dp4_groups
        #=while true
            dp4 = correct_permutation(comm)
            zero_sum_rows = vec(sum(dp4[1], dims=2) .== 0)
    
            # Skip if all rows would be filtered out
            if all(zero_sum_rows)
                continue
            end
    
            # Keep only non-zero rows
            dp4_filtered = dp4[1][.!zero_sum_rows, :]
    
            # Filter groups the same way (assuming groups is a vector)
            dp4_groups = groups[.!zero_sum_rows]
    
            # Check if we have at least 2 unique groups
            if length(unique(dp4_groups)) >= 2
                break
            end
        end=#

        #SIMPER analysis performed on each permutated matrix
        simp2 = simper(dp2[i], groups)
        simp3 = simper(dp3[i], groups)
        #simp4 = simper(dp4_filtered, dp4_groups)
        simp4 = simper(dp4[i], groups)

        #Storage of SIMPER results (ranked contribution to OAD) and conversion to percentage of SIMPER results
        sorted2 = sort(vec(simp2[1]), rev=true)
        sorted3 = sort(vec(simp3[1]), rev=true)
        sorted4 = sort(vec(simp4[1]), rev=true)
        
        sum2 = sum(sorted2)
        sum3 = sum(sorted3)
        sum4 = sum(sorted4)
        
        df2[i,:] = (sorted2/sum2) * 100
        df3[i,:] = (sorted3/sum3) * 100
        df4[i,:] = (sorted4/sum4) * 100

        
        #if isnan(sum(df4[i,:]))
        #    println("There are NaNs in the SIMPER results of the permuted matrix when the column sum is maintained, $i")
        #    println("$sorted4, $sum4, $(df4[i,:])")
        #end

    end

    dn2=hcat([sort(df2[:, j]) for j in 1:size(df2, 2)]...)
    dn3=hcat([sort(df3[:, j]) for j in 1:size(df3, 2)]...)
    dn4=hcat([sort(df4[:, j]) for j in 1:size(df4, 2)]...)
    #println("$dn4")

    up = floor(Int,0.975*Nperm) # ex for 100 permutations it will used 97
    lo = floor(Int,0.025*Nperm) # ex for 100 permutations it will used 2
    med = floor(Int,0.5*Nperm)

    ###The following is the calculation of E index

    # Ranked % of contribution to OAD of empirical and simulated profiles
    obs = Pourcent_Contribution
    #print("Empirical profile: $obs")
    Orange = dn4
    Blue = dn2
    Green = dn3

    VectorEcartCarreOrangeLog = zeros(Nperm)
    VectorEcartCarreGreenLog = zeros(Nperm)
    VectorEcartCarreBlueLog = zeros(Nperm)

    for i in 1:Nperm
        # Calculate the sum of square deviations for each profile
        # between the empirical profile and the simulated profiles
        SommeEcartCarreOrange = (Orange[i,:] .- obs).^2
        SommeEcartCarreGreen = (Green[i,:] .- obs).^2
        SommeEcartCarreBlue = (Blue[i,:] .- obs).^2
        
        # Log conversion of the sum of square deviations
    
        sum_orange = sum(SommeEcartCarreOrange)
        VectorEcartCarreOrangeLog[i] = log10(sum_orange + 1.0e-20)

        sum_green = sum(SommeEcartCarreGreen)
        VectorEcartCarreGreenLog[i] = log10(sum_green + 1.0e-20)

        sum_blue = sum(SommeEcartCarreBlue)
        VectorEcartCarreBlueLog[i] = log10(sum_blue + 1.0e-20)

        if VectorEcartCarreOrangeLog[i]==-Inf
            println("-Inf in Orange, $i")
        end
    end

    DataMeanCarreLog = DataFrames.DataFrame(
                        Orange = VectorEcartCarreOrangeLog,
                        Blue = VectorEcartCarreBlueLog,
                        Green = VectorEcartCarreGreenLog)
  
    final_result = Dict("EcartCarreLog" => DataMeanCarreLog,
                        "mat" => comm, 
                        "ContriPercentage" => Pourcent_Contribution,
                        "UpOrange" => dn4[up,:], 
                        "DownOrange" => dn4[lo,:], 
                        "MedOrange" => dn4[med,:],
                        "UpBlue" => dn2[up,:], 
                        "DownBlue" => dn2[lo,:], 
                        "MedBlue" => dn2[med,:],
                        "UpGreen" => dn3[lo,:], 
                        "DownGreen" => dn3[up,:], 
                        "MedGreen" => dn3[med,:],
                        "dnOrange" => dn4, 
                        "dnBlue" => dn2, 
                        "dnGreen" => dn3)
    return final_result
end

#A function to calculate DNCI for only two groups
#This function is a translation/adaptation of a function from the R package `DNCImper`, licensed under GPL-3.
#Original package and documentation available at: https://github.com/Corentin-Gibert-Paleontology/DNCImper
function DNCI_ses(comm::Matrix, groups::Vector, Nperm::Int=1000; count::Bool=false) #for presence-absence data only
    group=sort(unique(groups))
    
    # Check if the number of groups is equal to 2
    if length(group) != 2
        error("length(groups) must be 2")
    end

    results = PerSIMPER(comm, groups, Nperm; count)
    E = results["EcartCarreLog"]

    # Precompute statistics
    mean_blue = mean(E.Blue)
    std_blue = std(E.Blue, corrected=true)
    cv = abs(std_blue / mean_blue)

    # Handle cases when both dispersal and niche processes are so constraining that there's only one or very limited possible arrangement of species across sites when performing quasi-swap
    if mean_blue == -20 && std_blue == 0  # Case 1
        DNCI = NaN
        S_DNCI = NaN
        CI_DNCI = NaN
        status = "quasi_swap_permutation_not_possible"
    elseif std_blue == 0  # Case 2
        DNCI = NaN
        S_DNCI = NaN
        CI_DNCI = NaN
        status = "one_way_to_quasi_swap"
    elseif cv < 0.01  # Case 3
        DNCI = NaN
        S_DNCI = NaN
        CI_DNCI = NaN
        status = "inadequate_variation_quasi_swap"
    else  # Normal case
        # Regular DNCI calculation
        SES_d = (E.Orange .- mean_blue) ./ std_blue
        SES_n = (E.Green .- mean_blue) ./ std_blue

        DNCI = mean(SES_d) - mean(SES_n)
        S_DNCI = sqrt(std(SES_d, corrected=true)^2 + std(SES_n, corrected=true)^2)
        CI_DNCI = 2 * S_DNCI
        status = "normal"
    end                               
    
    metric = DataFrames.DataFrame(group1= group[1], 
    group2 = group[2], 
    DNCI = DNCI, 
    CI_DNCI = CI_DNCI, 
    S_DNCI = S_DNCI,
    status = status)

    return metric
end

end