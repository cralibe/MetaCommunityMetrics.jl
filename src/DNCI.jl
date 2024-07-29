using Pkg
Pkg.activate(".")

using Distributions
using Plots
using LinearAlgebra
using GaussianRandomFields
using Distances
using DataFrames
using CSV
using StatsBase
using DataStructures
using ProgressMeter
using Pipe: @pipe
using Combinatorics

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
# A function to create groups for each year
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

# A function to plot the groups
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


