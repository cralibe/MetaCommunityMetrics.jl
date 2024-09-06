using Pkg
Pkg.activate()

using CSV
using DataFrames
using Pipe

#Example Data
#Load sample data included in the package
sample_df = CSV.read(joinpath(@__DIR__, "..", "..", "data", "rodent_abundance_data.csv"), DataFrame, skipto=2)#read in the sample data
sample_df = stack(sample_df, Not(:Sampling_date_order, :Year, :Month, :Day, :plot), variable_name = :Species, value_name = :Abundance) #stack the dataframe


#A dataframe that idetifies the persisting species (The species that existed at least once in all time steps)
persist_sp_df=
@pipe sample_df |> 
groupby(_, :Species) |> #group the dataframe by every species
combine(_, :Abundance => sum => :Total_Abundance) |> #calculate the total abundance of every spesice within the whole sampling period
filter(row -> row[:Total_Abundance] > 0, _) #select rows that has a total abundance >0.

#A dataframe that idetifies non empty patch at every time step
non_empty_site = 
@pipe sample_df[in(persist_sp_df.Species).(sample_df.Species), :] |>
groupby(_, [:Sampling_date_order, :plot]) |>
combine(_, :Abundance => sum => :Total_Abundance) |>
filter(row -> row[:Total_Abundance] > 0, _)

#The data that contains the abundace and presence-absence of persisting species and non empty site at every time step
metacomm_df=
@pipe sample_df|>
_[in(persist_sp_df.Species).(_.Species), :]|> #keeping species with a total abundance >0 across all patches all time steps. 
innerjoin(_, non_empty_site, on = [:Sampling_date_order, :plot]) |> #keeping non empty site at every time step. 
select(_, Not(:Total_Abundance)) |> #remove the Total_Abundance column to avoid confusion
transform(_, :Abundance => ByRow(x->ifelse.(x> 0, 1, 0)) => :Presence) #create a presence-absence column from the abundance column

#Simulating coordinates for the sample data. This step is not nessary if: (1) you have the coordinates data, (2) you don't need to caluculating DNCI, or (3) you are doing the groupings for the DNCI by yourself.
#Set grid dimensions
num_rows = 4
num_columns = 6
num_plots = 24  # Total number of plots
# Define base coordinates (for example purposes, these are arbitrary)
base_latitude = 35.0
base_longitude = -110.0
# Define distance between plots in degrees (arbitrary small values for demonstration)
lat_spacing = 0.5
lon_spacing = 0.5
# Initialize DataFrame
patch_coord_df = DataFrame(plot = Int[], Latitude = Float64[], Longitude = Float64[])

# Generate grid coordinates
plot_count = 1
for i in 1:num_rows
    for j in 1:num_columns
        if plot_count > num_plots
            break
        end
        lat = base_latitude + (i - 1) * lat_spacing
        lon = base_longitude + (j - 1) * lon_spacing
        push!(patch_coord_df, (plot_count, lat, lon))
        plot_count += 1
    end
end

#Adding the coordinates to the metacomm_df
metacomm_df = @pipe metacomm_df |>
innerjoin(_, patch_coord_df, on = :plot) #joining the metacomm_df with the patch_coord_df

#Preparing the matric for calculating beta diveristy
#matrix with the species abundance
sample_matrix_abundance = @pipe metacomm_df |> 
select(_, Not(:Presence)) |>
filter(row -> row[:Sampling_date_order] == 50, _) |> # filter rows where Sampling_date_order == 50
unstack(_, :Species, :Abundance, fill=0) |> #convert it back to the wide format 
select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
Matrix(_) |> #convert the dataframe to a matrix
_[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero
#matrix with the species presence-absence data
sample_matrix_presence = @pipe metacomm_df |> 
select(_, Not(:Abundance)) |>
filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
Matrix(_) |> #convert the dataframe to a matrix      
_[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero

#Preparing the data for the DNCI analysis
total_presence_df=@pipe metacomm_df|>
groupby(_,[:Species,:Sampling_date_order])|>
combine(_,:Presence=>sum=>:Total_Presence)
#filter(row -> row[:Total_Presence] <= 1, _)|> 

#Remove singletons 
total_richness_df= @pipe metacomm_df|>
leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
filter(row -> row[:Total_Presence] > 1, _)|> 
groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
combine(_,:Presence=>sum=>:Total_Richness)

#Prepare a community matrix at t=50 for the DNCI analysis

comm= @pipe metacomm_df|>
leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
filter(row -> row[:Total_Presence] > 1, _) |>
filter(row -> row[:Sampling_date_order] == 50, _) |>
select(_, [:plot, :Species, :Presence]) |>
unstack(_, :Species, :Presence, fill=0) |>
select(_, Not(:plot)) |>
Matrix(_)

cluster_result = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness) 
        
presence_df = @pipe metacomm_df|>
leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
filter(row -> row[:Total_Presence] > 1, _) |>
filter(row -> row[:Sampling_date_order] == 50, _) |>
select(_, [:plot, :Species, :Presence]) |>
unstack(_, :Species, :Presence, fill=0) |>
leftjoin(_, cluster_result[50], on = [:plot => :Patch], makeunique = true) 