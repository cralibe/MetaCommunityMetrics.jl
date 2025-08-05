using Pkg
Pkg.activate(".")

using CSV
using DataFrames
using Pipe
using Random
using Distributions
using UnicodePlots

#Sample Data
#Load sample data included in the package
sample_df = CSV.read(joinpath(pkgdir(MetaCommunityMetrics), "data", "rodent_abundance_data.csv"), DataFrame, skipto=2)#read in the sample data
sample_df = stack(sample_df, Not(:Sampling_date_order, :Year, :Month, :Day, :plot), variable_name = :Species, value_name = :Abundance) #stack the dataframe


#A dataframe that idetifies the persisting species (The species that existed at least once in all time steps)
persist_sp_df=
@pipe sample_df |> 
groupby(_, :Species) |> #group the dataframe by every species
combine(_, :Abundance => sum => :Total_Abundance) |> #calculate the total abundance of every spesice within the whole sampling period
filter(row -> row[:Total_Abundance] > 0, _) #select rows that has a total abundance >0.

#A dataframe that idetifies non empty patch during the sampling period (The patches that has at least one species with positive total abundance)
non_empty_site = 
@pipe sample_df[in(persist_sp_df.Species).(sample_df.Species), :] |>
groupby(_, [:plot]) |>
combine(_, :Abundance => sum => :Total_Abundance) |>
filter(row -> row[:Total_Abundance] > 0, _)

#The data that contains the abundace and presence-absence of persisting species and non empty site at every time step
metacomm_df=
@pipe sample_df|>
filter(row -> row.Species in persist_sp_df.Species, _) |> # filter species with positive total abundance
innerjoin(_, non_empty_site, on = [:plot]) |> # join with non-empty sites
select(_, Not(:Total_Abundance)) |> # remove `Total_Abundance` to avoid confusion
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

#Generating environmental data
env_df = @pipe metacomm_df |>
groupby(_, [:Sampling_date_order, :plot]) |>
combine(x -> DataFrame(standardized_temperature = 0.0, standardized_precipitation = 0.0), _)


#Generate random environmental data
Random.seed!(123)
temp = rand(Normal(15.0, 2.0), 2808)
precip = rand(LogNormal(log(50.0) - 0.5*log(1 + (30.0/50.0)^2), sqrt(log(1 + (30.0/50.0)^2))), 2808)

#Transform precipitation to achieve normality
log_precip = log.(precip)

#Standardize both variables (z-score: subtract mean, divide by standard deviation)
temp_standardized = (temp .- mean(temp)) ./ std(temp)
precip_standardized  = (log_precip .- mean(log_precip)) ./ std(log_precip)

env_df.standardized_temperature .= temp_standardized
env_df.standardized_precipitation .= precip_standardized

#Checking the distribution of the environmental data
UnicodePlots.histogram(env_df.standardized_temperature, 
                         nbins = 10, 
                         title = "Temperature Histogram",  
                         width=20,  
                         height=5)

UnicodePlots.histogram(env_df.standardized_precipitation,
                         nbins = 10, 
                         title = "Precipitation Histogram",  
                         width=20,  
                         height=5)

#Adding the coordinates and environmental parameters to the metacomm_df
metacomm_df = @pipe metacomm_df |>
innerjoin(_, patch_coord_df, on = :plot) |>#joining the metacomm_df with the patch_coord_df
leftjoin(_, env_df, on = [:Sampling_date_order, :plot]) #joining the metacomm_df with the env_df

#Save the metacomm_df to a csv file
CSV.write(joinpath(pkgdir(MetaCommunityMetrics), "data", "metacomm_rodent_df.csv"), metacomm_df)


