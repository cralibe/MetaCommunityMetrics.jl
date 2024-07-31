# src/DNCI.jl

using ..Internal


"""
    create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int}) -> Dict{Int, DataFrame}

Create clusters for each unique time step in a dataset. Only presnece-absence data can be used.

# Arguments
- `time::Vector`: A vector indicating the time each sample was taken.
- `latitude::Vector`: A vector indicating the latitude of each sample.
- `longitude::Vector`: A vector indicating the longitude of each sample.
- `patch::Vector`: A vector indicating the spatial location (patch) of each sample.

# Returns
- `Dict{Int, DataFrame}`: A dictionary where each key is a unique time from the dataset and each value is a DataFrame for that time with an added `Group` column indicating the assigned cluster.

# Details
This function performs hierarchical clustering on geographical coordinates for each unique time step. It aims to balance the clusters based on the number of sites and species richness, ensuring that no group has less than five sites and that there are at least two groups. If conditions for clustering balance are not met (like groups having less than five sites or only one group), it iteratively adjusts the clusters by reassigning sites to improve group balance.

# Example
```julia
using CSV, DataFrames
using DataFramesMeta
using Pipe: @pipe

sample_df = @pipe CSV.read("data/rodent_abundance_data.csv", DataFrame; header=true) |>#read in the sample data
            select(_, Not(:Column1))|> #select the columns 
            stack(_, Not(:Sampling_date_order, :Year, :Month, :Day, :plot), variable_name = :species, value_name = :abundance)

preped_data= @pipe sample_df |>         
@transform(_, :presence = ifelse.(:abundance .>= 1, 1, 0)) |>
groupby(_, [:Sampling_date_order, :species]) |> 
combine(_,:presence=>sum=>:total_presence) |>
filter(row -> row[:total_presence] !<= 1, _) |> #remove singletons (species occurring at one site only)
leftjoin(_,sample_df, on = [:Sampling_date_order, :species]) |> #join the data back to the original data


# Generate random latitude and longitude values
n = nrow(sample_df)
latitude = rand(35.0:0.01:36.0, n) # Adjust range as needed
longitude = rand(-120.0:0.01:-119.0, n) # Adjust range as needed

# Add the coordinates to the DataFrame
sample_df[:, :latitude] = latitude
sample_df[:, :longitude] = longitude

# Create groups for each time step
grouped_data = create_clusters(sample_df.Sampling_date_order, sample_df.latitude, sample_df.longitude, sample_df.plot)

```
"""
function create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int})
    grouping_dict = Dict{Int, DataFrame}()

    #Create a DataFrame
    df = DataFrame(Time = time, Latitude = latitude, Longitude = longitude, Patch = patch)

    for t in unique(df[:, :Time])
        subset_df = filter(row -> row[:Time] == t, df)
        coordinates = select(subset_df, [:Latitude, :Longitude])
        distances = Distances.pairwise(Euclidean(), Matrix(coordinates), dims=1)
        
        num_clusters = max(div(length(unique(subset_df.Patch)), 5), 2) # Ensure that the number of clusters is at least 2 and each group has at least 5 sites
        condition_met = false

        while !condition_met && num_clusters > 1
            agglo_result = hclust(distances, linkage=:complete)
            assignments = cutree(agglo_result, k=num_clusters)
            subset_df.Group = assignments
            subset_df = check_condition_and_fix(subset_df)
            condition_met = check_conditions(subset_df)
            
            if !condition_met
                num_clusters -= 1
                if num_clusters < 2
                    println("Warning: Cluster count fell below 2, which is not permissible for clustering. Groups assigned as missing.")
                    subset_df.Group .= missing
                end
            end
        end

        grouping_dict[t] = subset_df
    end
    if length(grouping_dict) != length(unique(df.Time))
        println("Warning: Some time steps are missing!!!")
    end    
    return grouping_dict
end

"""
    plot_clusters(grouped_data::DataFrame)

This function plots the clustering result from the function create_clusters() on a scatter plot according to geographic coordinates.
Each point represents a location, and points are colored based on their group/cluster ID.

# Arguments
- `grouped_data::DataFrame`: A DataFrame containing the columns `Latitude`, `Longitude`, and `Group`. 
    - `Latitude`: A column of latitude coordinates.
    - `Longitude`: A column of longitude coordinates.
    - `Group`: A column indicating the cluster ID for each point.

# Returns
Returns
- This function displays a scatter plot of the clusters.

# Details
- The function extracts the `Latitude`, `Longitude`, and `Group` columns from the provided DataFrame.
- Unique cluster IDs are assigned numeric identifiers.
- A color palette is generated using distinguishable colors to differentiate the clusters.
- The scatter plot is created with points colored by their cluster ID.

# Example
```julia
using DataFrames
# Create a sample DataFrame
data = DataFrame(Latitude = [34.05, 36.16, 40.71, 34.05],
                 Longitude = [-118.24, -115.15, -74.01, -118.24],
                 Group = ["A", "B", "A", "B"])

# Plot clusters
plot_clusters(data)
```
"""
function plot_clusters(grouped_data::DataFrame)
    # Extract latitude, longitude, and cluster columns
    latitudes = grouped_data.Latitude
    longitudes = grouped_data.Longitude
    cluster_ids = grouped_data.Group

    # Get unique cluster IDs and assign numeric identifiers
    unique_clusters = unique(cluster_ids)
    cluster_map = Dict(cluster => i for (i, cluster) in enumerate(unique_clusters))

    # Assign numeric identifiers to cluster IDs
    numeric_ids = [cluster_map[cluster] for cluster in cluster_ids]

    # Define color palette for clusters (you can modify this or use any other color palette)
    colors = distinguishable_colors(length(unique_clusters), colorant"blue")

    # Plot the points, color by cluster ID
    scatter(longitudes, latitudes, marker_z=numeric_ids,
            xlabel="Longitude", ylabel="Latitude", title="Cluster Visualization $(subset_df.Time[1])",
            legend=false, markerstrokecolor=:black, markerstrokewidth=0.5,
            markersize=5, color=colors[numeric_ids], label=false)
end

"""
    DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000, count::Bool=true)

Calculates the dispersal-niche continuum index (DNCI) for multiple groups, a metric proposed by Vilmi(2021) (doi: 10.1111/ecog.05356).

# Arguments
- `comm::Matrix`: A presence-absence data matrix where rows represent observations (e.g., sites or samples) and columns represent species.
- `groups::Vector`: A vector indicating the group membership for each row in the `comm` matrix.
- `Nperm::Int=1000`: The number of permutations for significance testing. Default is 1000.
- `count::Bool=true`: A flag indicating whether the numeber of permutations is printed. Default is `false`.

# Returns
- `DataFrame`: A DataFrame containing the DNCI results for each pair of groups.


# Examples
```julia
# Example usage of DNCI_multigroup
comm = [1 0 0 1 0;
        1 1 0 0 0;
        0 1 1 0 0;
        0 0 1 1 1;
        1 0 0 0 1;
        0 1 1 0 1]

groups = ["A", "A", "B", "B", "C", "C"]
Nperm = 1000
count = true

result = DNCI_multigroup(comm, groups, Nperm, count)
println(result)

"""
function DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000, count::Bool=true) #for presence-absence data only
    group_combinations = collect(combinations(unique(sort(groups)),2))

    ddelta = DataFrame()

    for i in 1:size(group_combinations,1)
        # Create an empty dictionary to hold the split data
        splitx = Dict()

        # Assume comm is a matrix and groups is a vector indicating the group for each row in comm
        for row in 1:size(comm, 1)  # Iterate over the rows of the matrix
            group = groups[row]  # Identify the group of the current row
            current_row = comm[row, :]   # Extract the entire row as a vector, make it a 1-row matrix
            row_matrix = reshape(current_row, 1, length(current_row))  # Reshape row to be a 1-row matrix

            # Check if the group already exists in the dictionary
            if haskey(splitx, group)
            # Vertically concatenate the new row matrix with the existing matrix for this group
            splitx[group] = vcat(splitx[group], row_matrix)
            else
            # Initialize the group with the current row as a 1-row matrix
                splitx[group] = row_matrix
            end
        end

        # Safely access and concatenate matrices from two groups
        if haskey(splitx, group_combinations[i][1]) && haskey(splitx, group_combinations[i][2])
            paired_x = vcat(splitx[group_combinations[i][1]], splitx[group_combinations[i][2]])
        else
            error("One of the groups does not exist in the dictionary.")
        end
        # Calculate a logical array indicating non-zero sum columns
        non_zero_sum_columns = sum(paired_x, dims=1) .!= 0
        # Convert logical index to actual column indices
        column_indices = findall(x -> x, non_zero_sum_columns[:]) 
    
        # Subset the matrix using these indices
        paired_x = paired_x[:,column_indices]

        group_pair = vcat(
        fill(group_combinations[i][1], size(splitx[group_combinations[i][1]], 1)),
        fill(group_combinations[i][2], size(splitx[group_combinations[i][2]], 1)))

       

        DNCI_result = DNCI_ses(paired_x, group_pair, Nperm, count)

        append!(ddelta, DNCI_result, promote=true)
    end
    return ddelta
end
  