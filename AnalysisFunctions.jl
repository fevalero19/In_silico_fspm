using CSV
using DataFrames
using Plots
using StatsBase
import Base.Filesystem.readdir
using HypothesisTests
using GLM
using LinearAlgebra


using Distances
using Clustering
using StatsPlots

using Statistics

"""
NB: SENSITIVITY ANALYSIS MUST BE DONE PRIOR TO THESE FUNCTIONS
    - THE PATH SHOULD CONTAIN SA_output/outputs/"

This script is meant to be used for parts of the calibration where the sensitivity of the parameters
    is again of interest. For instance when the precision of the estimated parameters is evaluated.

Most functions are the same as in AnalyzeSensitivities2. However some have been ##ADAPTED
"""



function retrieve_parameter(k::String; plus::Bool=true, spec = 1)
    v = p[spec][k]
    if isa(v, Int64)
        if plus
            v = Int64(ceil(v*1.1))
        else
            v = Int64(floor(v*0.9))
        end
        # Did not work for nleaflets.
    elseif isa(v, Float64) 
        if plus
            v = v*1.1
        else
            v = v*0.9
        end
    else
        println("$k")
    end
    return v
end

function obtain_parameter_names()
    file_names = readdir("data/")
    # Create empty vectors to store entries containing "D1" and "S2"
    p1 = String[]
    p2 = String[]

    # Filter and populate p1 and p2
    for entry in file_names
        if contains(entry, "D1") && contains(entry, "field") && !contains(entry, "Control")
            new_entry = entry[12:end-7]
            if contains(new_entry, "+D")
                error("A double file has been observed, Please fix file : $entry")
                # This line was added for the potential_biomass_flower file. It seems that this simulation was run twice, probably because of the moment that Random.seed was found active.
                    # The potential_biomass_flower plants2 occured 2 times. The newest was chosen.
                # solution: the older file was removed

                # Very weird
                    # In data_plants3half_max_sheat+S2.txt there are 6 wheat plants missing at timestep 48. Just that, the first 6 triticale plants are missing
                    # This created a weird bug where all of a sudden the time vector was converted into a PooledArrays.PooledVector{String, UInt32, Vector(Uint32)}
                # solution: the data_plants2half_max_sheath+S2.txt was copied and renamed.
                    # CHECK IMMEDIATLY WHETHER half_max_sheath became important. 
            end
            push!(p1, new_entry)
        elseif contains(entry, "S2") && contains(entry, "field") && !contains(entry, "Control")
            new_entry = entry[12:end-7]
            
            # The field3endosperm_release_speed-S2 was not yet done, which was why the field2-S2 was copied and renamed as -S2.
                # Same for plants
            
            if contains(new_entry, "+S")
                error("A double file has been observed, Please fix file : $entry")
            end
            push!(p2, new_entry)
        end
    end

    p1_u = unique(p1)
    p2_u = unique(p2)
    must_filter = false
    for n in p1_u
        if count(x->x==n, p1) < 6
            println("The parameter $n D1 has less than 6 simulations")
            must_filter = true
        end
    end
    if must_filter
        p1_u = [n for n in p1_u if count(x -> x == n, p1) >= 6]
        must_filter = false
    end
    for n in p2_u
        if count(x->x==n, p2) < 6
            println("The parameter $n S2 has less than 6 simulations")
            must_filter = true
        end
    end 
    if must_filter  # Solution for now is to remove the parameter max_branch_angle from species 2 analysis. Until the simulations are complete.
            # If then the number is still wrong, a new set of simulations will have to be done
        p2_u = [n for n in p2_u if count(x -> x == n, p2) >= 6]
        must_filter = false
    end



    return p1_u, p2_u
end


## ADAPTED!!
    # This method is much faster (12.7x) than the original function. 
function obtain_output_locations()
    y_field = unique([p[7:end-7] for p in filter(x-> contains(x,"field"), readdir("SA_output/outputs"))])
    y_plants = unique([p[7:end-7] for p in filter(x-> contains(x,"plant"), readdir("SA_output/outputs"))])
    return y_plants, y_field
end

function obtain_control_files(level)
    if level == "field"
        filenames = ["data_fieldControlD1_1_1",
                     "data_fieldControlD1_2_1",
                     "data_fieldControlD1_3_1",
                     "data_fieldControlD2_1",
                     "data_fieldControlD2_2",
                     "data_fieldControlD2_3"]
    elseif level == "plant"
        filenames = ["data_plantsControlD1_1_1",
                    "data_plantsControlD1_2_1",
                    "data_plantsControlD1_3_1",
                    "data_plantsControlD2_1",
                    "data_plantsControlD2_2",
                    "data_plantsControlD2_3"]
    else
        error("$level is not a correct level of output")
    end
    return filenames
end

function output_level(output)
    y_p, y_f = obtain_output_locations()
    if output in y_p
        return "plant"
    elseif output in y_f
        return "field"
    else
        error("Output $output could not be found in plant or field file")
    end
end

function get_control_data(output, species; Stats = "mean")
    level = output_level(output)
    if level == "field"
        DF = CSV.read("data/data_fieldControlD1_1.txt", DataFrame)
        timeseries = unique(DF.time)
        species_str =string(species)
        colnames = ["$output-$species_str-t$t" for t in timeseries]
        row = [filter(row -> row.time == t && row.species == species, DF)[!,output] for t in timeseries]
        controlDF = DataFrame(row, colnames)
        for c in obtain_control_files(level)[2:end] # The first txt file had all ready been read
            DF = CSV.read("data/$c.txt", DataFrame)
            row = [filter(row -> row.time == t && row.species == species, DF)[!,output] for t in timeseries]
            push!(controlDF, hcat(row...))
        end
    elseif level == "plant"
        DF = CSV.read("data/data_plantsControlD1_1_1.txt", DataFrame; delim = "\t")
        timeseries = unique(DF.time)
        species_str =string(species)
        colnames = ["$output-$species_str-t$t" for t in timeseries]
        if Stats == "mean"
            row = [[mean(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
        elseif Stats == "std"
            row = [[std(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
        else
            println("$Stats not used, provided the mean")
            row = [[mean(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
        end
        controlDF = DataFrame(row, colnames)
        for c in obtain_control_files(level)[2:end] # The first txt file had all ready been read
            DF = CSV.read("data/$c.txt", DataFrame; delim = "\t")
            if Stats == "mean"
                row = [[mean(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
            elseif Stats == "std"
                row = [[std(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
            else
                println("$Stats not used, provided the mean")
                row = [[mean(filter(row -> row.time == t && row.species == species, DF)[!,output])] for t in timeseries]
            end
            push!(controlDF, hcat(row...))
        end
    end  
    return controlDF
end

function get_sensitivity_data(parameter, speciesP, output, speciesY; Stats = "mean")
    global p
    level = output_level(output)
    if speciesP == 1
        s = "D1"
    elseif speciesP == 2
        s = "S2"
    else
        error("$speciesP is not a correct value for species")
    end

    if level == "field"
        data_for_sensitivity = get_control_data(output, speciesY)
        timeseries = 1:size(data_for_sensitivity)[2]
        # Obtain the species present in the control file
        insertcols!(data_for_sensitivity, :parameter_value => repeat([p[speciesP][parameter]],size(data_for_sensitivity)[1]))
    
        for k in 1:3
            Field_data = CSV.read("data/data_field$k$parameter+$s.txt", DataFrame)
            y = [filter(row -> row.time == t && row.species == speciesY, Field_data)[!,output] for t in timeseries]
            row = hcat(y..., retrieve_parameter(parameter; plus = true, spec = speciesP))
            push!(data_for_sensitivity, row)
        end
        for k in 1:3
            Field_data = CSV.read("data/data_field$k$parameter-$s.txt", DataFrame)
            y = [filter(row -> row.time == t && row.species == speciesY, Field_data)[!,output] for t in timeseries]
            row = hcat(y..., retrieve_parameter(parameter; plus = false, spec = speciesP))
            push!(data_for_sensitivity, row)
        end
        return data_for_sensitivity
    elseif level == "plant"
        data_for_sensitivity = get_control_data(output, speciesY)
        timeseries = 0:size(data_for_sensitivity)[2]-1
        # Obtain the species present in the control file
        insertcols!(data_for_sensitivity, :parameter_value => repeat([p[speciesP][parameter]],size(data_for_sensitivity)[1]))
        for k in 1:3
            #println("data/data_plants$k$parameter+$s.txt")
            DF = CSV.read("data/data_plants$k$parameter+$s.txt", DataFrame; delim = "\t")
            if Stats == "mean"
                y = [[mean(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            elseif Stats == "std"
                y = [[std(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            else
                println("$Stats not used, provided the mean")
                y = [[mean(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            end
            row = hcat(y..., retrieve_parameter(parameter; plus = true, spec = speciesP))
            push!(data_for_sensitivity, hcat(row...))
        end
        for k in 1:3
            DF = CSV.read("data/data_plants$k$parameter-$s.txt", DataFrame; delim = "\t")
            if Stats == "mean"
                y = [[mean(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            elseif Stats == "std"
                y = [[std(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            else
                println("$Stats not used, provided the mean")
                y = [[mean(filter(row -> row.time == t && row.species == speciesY, DF)[!,output])] for t in timeseries]
            end
            row = hcat(y..., retrieve_parameter(parameter; plus = false, spec = speciesP))
            push!(data_for_sensitivity, hcat(row...))
        end
        return data_for_sensitivity
    else
        error("Not correct level found $level")
    end
    
    
    return 
end

function calculateSensitivity(parameter_names, speciesP, output, speciesY)
    S = zeros(length(1:115), length(parameter_names))
    for (i,p) in enumerate(parameter_names)    # Go through parameters
        SAdf = get_sensitivity_data(p, speciesP, output, speciesY)
        uppervalues = filter(row -> row.parameter_value == maximum(SAdf.parameter_value), SAdf)
        lowervalues = filter(row -> row.parameter_value == minimum(SAdf.parameter_value), SAdf)
            # Test if + and - are significantly different
        for j in 1:115
            result = HypothesisTests.EqualVarianceTTest(uppervalues[!,j], lowervalues[!, j])
            if pvalue(result) < 0.05    # We can assume that the parameter had a significant influence on the outcome
                S4 = coef(GLM.lm(@formula(y~x), DataFrame(x=SAdf[!,"parameter_value"], y=SAdf[!,j])))[2]
                S[j,i] = S4
            end
        end
    end
    S = DataFrame(S, parameter_names)
    return S
end

function normalizeS(S::DataFrame, speciesP, output::String, speciesY)
    ps = names(S)
    avg_p = [p[speciesP][k] for k in ps]
    control_df = get_control_data(output, speciesY)
    avg_y = [mean(c) for c in eachcol(control_df)]
    S = S |>Matrix   # Convert S to a Matrix
    NS = S .* (1 ./avg_y * avg_p')
    NS = replace(NS, NaN => 0)
    NS = DataFrame(NS, ps)
    return NS
end


function construct_corrMatrix(mJ)
    """
    This function is very fragile to problems in mJ
    """
    cov = ones(size(mJ)[2],size(mJ)[2])
    plot_heat = false
    try
        cov = LinearAlgebra.inv(mJ'*mJ)
        plot_heat = true
        if !all(i->i>=0,diag(cov))
            cov = ones(size(mJ)[2],size(mJ)[2])
            plot_heat = false
        end
    catch
       println("Correlation matrix could not be setup") 
    end
    corr = deepcopy(cov)
    std_deviations = sqrt.(diag(cov))
    number_parameters = length(diag(cov)) #!!
    for i in 1:number_parameters, j in 1:number_parameters
        corr[i, j] /= (std_deviations[i] * std_deviations[j])
    end
    if plot_heat
        display(StatsPlots.heatmap(abs.(corr)))
    end
    return corr
end


function cluster_param(corr)
    """
    Returns for each parameter the index of the cluster they are in.
    """
    if !issymmetric(corr)   # There are numerical errors in calculating the correlation matrix
        dgts = 1
        while issymmetric(round.(corr, digits = dgts))  # This rounds the correlation matrix at later digits until it is no longer symmetric
            dgts += 1
        end
        if dgts<3
            println("The non-symmetry of the correlation matrix can not be explained by numerical error. The number of digits was: $dgts")
        end
        corr = round.(corr, digits = dgts-1)    # Then the number 
    end
    dist = 1 ./corr
    hc = Clustering.hclust(abs.(dist), linkage=:single) # Single linkage to ensure that the closest parameter correlation lies at the cutoff height
    a = StatsPlots.plot(hc) # Plot dendogram
    threshold = 1/critical_corr # Adapt to how the distance is calculated
                        # When multiple outputs are considered, this threshold can be 0.75
                        # When only 1 output is considered, this threshold should be at most 0.65, but best 0.49.
                        # The identifiability test quantifies the structural identifiability
                        # Critical corr
    hline!([threshold], color =:red, label = "threshold")
    display(a)
    clusters = cutree(hc, h = threshold)
    return clusters
end

function param_calib(NmJacB, clusters)
    # Example data
    param_names = copy(parameter_names)
    jacobian = NmJacB  # Example Jacobian matrix
    cluster_vector = clusters # Cluster assignments for each parameter
    # Find unique clusters
    unique_clusters = unique(cluster_vector)

    # Initialize dictionaries to store the parameter names and values for the highest in each cluster
    highest_values = Dict{Int, Float64}()
    highest_param_names = Dict{Int, String}()

    # Iterate through clusters
    for cluster in unique_clusters
        # Find indices of parameters in the cluster
        parameter_indices = findall(p -> p == cluster, cluster_vector)
        if length(parameter_indices)==1
            parameter_index = parameter_indices[1]
            max_value = maximum(abs.(jacobian[:,parameter_index])) # To indicate that it was alone
            if max_value > eta1
                highest_param_names[cluster] = param_names[parameter_index]
                highest_values[cluster] = max_value
            end
        else
            # Extract the values from the Jacobian for parameters in the cluster
            cluster_values = jacobian[:, parameter_indices]
            
            # Find the parameter index with the maximum value in the cluster
            max_value, max_index = findmax(abs.(cluster_values))
            # Get the corresponding parameter name
            if max_value > eta1
                param_name = param_names[parameter_indices[max_index[2]]] # The max_index is a Cartesian that includes the moment of the maximum, only the second indicates the parameter
            
                # Store the parameter name and its highest value in the dictionary
                highest_param_names[cluster] = param_name
                highest_values[cluster] = max_value
            end
            
        end
    end

    # Determine the names of the other parameters
    calibratable_parameters = [p for p in values(highest_param_names)]

    other_param_names = setdiff(param_names, calibratable_parameters)

    # Print the results
    println("Parameters with the highest value in each cluster:")
    for (pN,pV) in zip(values(highest_param_names), values(highest_values))
        println("$pN (Value: $pV)")
    end

    println("\nOther parameters:")
    println(other_param_names)


    #setdiff(param_names, other_param_names)
    return calibratable_parameters, values(highest_values), other_param_names
end

function read_NS(output, speciesY, parameter::Vector{String}, speciesP::Int)
    level = output_level(output)
    if level == "plant"
        NS = CSV.read("SA_output/outputs/plant_$output$speciesY-$speciesP.csv", DataFrame)
    elseif level == "field"
        NS = CSV.read("SA_output/outputs/field_$output$speciesY-$speciesP.csv")
    else
        error("The level of output $output could not be found. The output is likely not to exist in the files of SA_output/outputs")
    end
    NS = NS[!,parameter]
    return NS
end

foo = read_NS("height", 2, ["endosperm", "plastochron", "specific_internode_length", "potential_biomass_internode"], 2)