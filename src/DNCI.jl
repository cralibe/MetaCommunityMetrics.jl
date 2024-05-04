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


###Dispersal niche continnum index
#load in all the functions related to the dispersal niche continuum index
include("DNCI.jl")
include("Clustering.jl")

function DNCI_with_clustering(presence_df, patch_coord_df) #the patch_coord_df need to specific a column 
    #with "North" and "South".
    ##Preparing the data for clustering
    #create a vector to contains singletons (species occurring at one site only)
    singletons_df=@pipe presence_df|>
    groupby(_,[:Species,:Time])|>
    combine(_,:Presence=>sum=>:Total_Presence)|>
    filter(row -> row[:Total_Presence] == 1, _)|>
    select(_,:Species)|>
    unique(_)|>
    Matrix(_)

    #Remove singletons and left join with the XY coordinates of the patches and the two ecoregions
    presence_xy_df=
    @pipe presence_df|>
    filter(row -> !(row[:Species] in singletons_df), _)|> #remove singletons
    leftjoin(_,patch_coord_df,on=:Patch)|>
    groupby(_,[:Patch,:Time,:Longitude,:Latitude, :LN])|>
    combine(_,:Presence=>sum=>:Total_Richness)

    #Seperate the data according to the two ecoregions
    presence_xy_df_North = presence_xy_df[presence_xy_df.LN .== "North", :] 
    presence_xy_df_South = presence_xy_df[presence_xy_df.LN .== "South", :]

    ##Clustering
    northern_groups_dict = create_clusters(presence_xy_df_North)
    southern_groups_dict = create_clusters(presence_xy_df_South)
    all_patches_groups_dict = create_clusters(presence_xy_df)

    #Concatenate the groupings dictionaries into dataframes
    #For north and south, respectively
    northern_groups_df=vcat(values(northern_groups_dict)...)
    southern_groups_df=vcat(values(southern_groups_dict)...)
    two_ecoregions_groups_df = vcat(northern_groups_df,southern_groups_df) #grouping based on the two ecoregions
    two_ecoregions_groups_df.Group[two_ecoregions_groups_df.LN .== "South"] .+= maximum(unique(northern_groups_df.Group)) #to avoid overlapping group numbers between the two ecoregions.
    #For all patches
    all_patches_groups_df=vcat(values(all_patches_groups_dict)...)

    #Checking if the two ecoregions sharing a similar species pool.
    #=two_group=@pipe presence_df|>
    filter(row -> !(row[:Species] in singletons_df), _)|> #remove singletons
    leftjoin(_,patch_coord_df,on=:Patch)|>
    filter(row -> row[:Presence] >0, _)|>#remove species that are not present in the sites
    groupby(_,:LN)

    set1=Set(two_group[1].Species)
    set2=Set(two_group[2].Species)
    shared_elements = intersect(set1, set2) #find the shared species between the two ecoregions
    #they shared 93 species in the observed data.=#

    ##Preparing data for the DNCI calculation
    #For the north and south
    two_ecoregions_groups_presence_df=@pipe presence_df |>
    unstack(_, :Species, :Presence)|> #unstack the data to create a presence-absence matrix
    transform!(_, names(_) .=> (c -> coalesce.(c, 0)) .=> names(_))|> #replace missing values with zeros
    leftjoin(_,two_ecoregions_groups_df,on= [:Patch, :Time])|>#left join the data with the grouping data
    select(_, :Patch, :Time, :Group, :LN, Not([:Patch, :Time, :Group, :LN, :Longitude, :Latitude, :Total_Richness]))
    #For all patches
    all_patches_groups_presence_df=@pipe presence_df |>
    unstack(_, :Species, :Presence)|> #unstack the data to create a presence-absence matrix
    transform!(_, names(_) .=> (c -> coalesce.(c, 0)) .=> names(_))|> #replace missing values with zeros
    leftjoin(_,all_patches_groups_df,on= [:Patch, :Time])|>#left join the data with the all patches grouping data
    select(_, :Patch, :Time, :Group, :LN, Not([:Patch, :Time, :Group, :LN, :Longitude, :Latitude, :Total_Richness]))


    #Calculate the DNCI
    #North
    DNCI_result_north = Dict{Int, DataFrame}()
    @showprogress 1 "Calculating DNCI for patches in the northern region..." for t in unique(two_ecoregions_groups_presence_df.Time)
        println("Calculating DNCI for time $t")
            df = @pipe two_ecoregions_groups_presence_df|>
            filter(row -> row[:LN] == "North", _)|>
            filter(row -> row[:Time] == t, _)
            
            comm = Matrix(df[:,6:end])
            groups = df.Group

        if size(comm, 1) != length(groups)
            error("Error: The number of rows in the data frame does not match the number of groups at time $t")
        end
            
            DNCI_result_north[t] = DNCI_multigroup(comm,groups,t)
    end
    #South
    DNCI_result_south = Dict{Int, DataFrame}()
    @showprogress 1 "Calculating DNCI for patches in the southern region..." for t in unique(two_ecoregions_groups_presence_df.Time)
        println("Calculating DNCI for time $t")
            df = @pipe two_ecoregions_groups_presence_df|>
            filter(row -> row[:LN] == "South", _)|>
            filter(row -> row[:Time] == t, _)
            
            comm = Matrix(df[:,6:end])
            groups = df.Group

        if size(comm, 1) != length(groups)
            error("Error: The number of rows in the data frame does not match the number of groups at time $t")
        end
            
            DNCI_result_south[t] = DNCI_multigroup(comm,groups,t)
    end
    #All patches
    DNCI_result_all_patches = Dict{Int, DataFrame}()
    @showprogress 1 "Calculating DNCI for all patches..." for t in unique(all_patches_groups_presence_df.Time)
        println("Calculating DNCI for time $t")
        df = @pipe all_patches_groups_presence_df|>
            filter(row -> row[:Time] == t, _)

            
        comm = Matrix(df[:,6:end])
        groups = df.Group

        if size(comm, 1) != length(groups)
            error("Error: The number of rows in the data frame does not match the number of groups at time $t")
        end

        DNCI_result_all_patches[t] = DNCI_multigroup(comm,groups,t)
    end

    #Concatenate the final result into dataframes
    DNCI_result_north_df = vcat(values(DNCI_result_north)...)
    DNCI_result_north_df.region = fill("North", size(DNCI_result_north_df, 1))
    DNCI_result_south_df = vcat(values(DNCI_result_south)...)
    DNCI_result_south_df.region = fill("South", size(DNCI_result_south_df, 1))
    DNCI_result_all_patche_df = vcat(values(DNCI_result_all_patches)...)
    DNCI_result_all_patche_df.region = fill("All", size(DNCI_result_all_patche_df, 1))

    DNCI_final_result= vcat(DNCI_result_north_df, DNCI_result_south_df, DNCI_result_all_patche_df)

    #Output the DNCI summary
    DNCI_summary=DataFrames.DataFrame(mean_DNCI_north = mean(DNCI_result_north_df.DNCI),
    min_DNCI_north = minimum(DNCI_result_north_df.DNCI),
    max_DNCI_north = maximum(DNCI_result_north_df.DNCI),
    mean_DNCI_south = mean(DNCI_result_south_df.DNCI),
    min_DNCI_south = minimum(DNCI_result_south_df.DNCI),
    max_DNCI_south = maximum(DNCI_result_south_df.DNCI),
    mean_DNCI_all = mean(DNCI_result_all_patche_df.DNCI),
    min_DNCI_all = minimum(DNCI_result_all_patche_df.DNCI),
    max_DNCI_all = maximum(DNCI_result_all_patche_df.DNCI))


   return DNCI_summary, DNCI_final_result
end




#pairwise_df=load("/home/jenny/phyto/niche_overlap_index_output/niche_overlap_matrix.jld2","pairwise_df")
