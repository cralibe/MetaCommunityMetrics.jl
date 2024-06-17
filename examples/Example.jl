#### A series of examples for using the functions provided by the MetaCommunity.jl
###Loading the required packages
using Pkg
Pkg.activate(.)
using CSV, DataFrames
using Pipe: @pipe
using MetaCommunity

##Suggested data wrangling before using any functions from the MetaCommunity package
#Load sample data included in the package and calculate beta diversity
sample_df = CSV.read("data/rodent_abundance_data.csv", DataFrame, skipto=2)#read in the sample data
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

first(metacomm_df, 5)#View the first five rows of the metacomm_df
#Simulating coordinates for the sample data. This step is not nessary if: (1) you have the coordinates data, (2) you don't need to caluculating DNCI, or (3) you are doing the groupings for the DNCI by yourself.
# Set grid dimensions
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

# Display the Coordinates DataFrame
println(patch_coord_df)

##Example 1: Calculating Beta Diversity at one time step (either using the abundance data or the presence-absence data)

#Preparing the matric for calculating beta diveristy
#matrix with the species abundance
sample_matrix_abundance = @pipe metacomm_df |> 
    select(_, Not(:Presence)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |> # filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Abundance, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero
#matrix with the species presence-absence data
sample_matrix_presence = @pipe metacomm_df |> 
    select(_, Not(:Abundance)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix      
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero
  
#calculate beta diversity for abundance data, set quant=true 
quant_result = beta_diversity(sample_matrix_abundance; quant=true) 

##Calculate beta diversity for presence-absence data, set quant=false
#Since this function can convert the data to binary data before calculating the beta diversity if abundance data is provided
#There are two methods to calculate beta diversity for presence-absence data

#First Method, using abundance data
binary_result_1 = beta_diversity(sample_matrix_abundance; quant=false)

#Second Method, using presence-absence data
binary_result_2 = beta_diversity(sample_matrix_presence; quant=false)

##Example 2: Calculating the mean spatial beta diversity (averaged across all time point) of a metacommunity
#We will use the metacomm_df here
first(metacomm_df, 5)#View the first five rows of the metacomm_df
mean_spatial_beta_div_result = mean_spatial_beta_div(metacomm_df.Abundance, metacomm_df.Sampling_date_order, metacomm_df.plot, metacomm_df.Species; quant=true)


##Example 3: Calculating the mean temporal beta diveristy (average across all patches) if a metacommunity

##Example 4: Calculating the variability metrics of a metacommunity

##Example 5: Calculating the occpied patches proportion

##Example 6: Calculating the niche overlap index

##Example 7: Calculating the Dispersal_Niche Contiuum Index (DNCI)

df = DataFrames.DataFrame(
        N=N,
        Time=Time,
        Patch=Patch,
        Species=Species)


test = @pipe metacomm_df |> 
    select(_, Not(:Presence)) |>
    filter(row -> row[:Year] == 2013, _) |>
    filter(row -> row[:Species] == "BA", _)

sum(test.Abundance)