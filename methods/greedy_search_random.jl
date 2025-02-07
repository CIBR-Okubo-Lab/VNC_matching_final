include("./MatchingFunctions.jl")
using ArgParse
using CSV
using DataFrames
using DataFramesMeta
using Dates
using JLD2
using Logging
using Random
using Statistics
using Wandb
using .MatchingFunctions

function main()
    # get time to synchronize solutions and Wandb
    now_time = string(Dates.format(now(), "yyyymmdd_HHMMSS"))
    parser = ArgParseSettings(description="The settings of Random Swap female edge algorithm")
    @add_arg_table! parser begin
        "--init_sol"
        help = "name of the initial solution file (without extension)"
        arg_type = String
        default = "5765447_20250107_101550.csv"

        "--seed"
        help = "random seed"
        arg_type = Int
        default = 392

        "--max_iter"
        help = "number of iterations to run"
        arg_type = Int
        default = 10000000
    end
    # init the random seed and file path
    args = parse_args(parser)
    Random.seed!(args["seed"])
    path_male = joinpath("..", "data", "male_connectome_graph.csv")
    path_female = joinpath("..", "data", "female_connectome_graph.csv")
    path_sol = joinpath("..", "sol", args["init_sol"])
    path_sub_graph = "./subgraphs.jld2"
    # read the male, female connection data and initial solution file
    df_male, df_female, df_sol = read_data(path_male, path_female, path_sol)
    # get the list of unique male and female nodes
    male_nodes = df_sol."Male Node ID"
    female_nodes = df_sol."Female Node ID"
    """read or get the subgraphs file,
       which is a dictionary of "Male Node"-> sub_graph DataFrame
    """
    subgraphs = Dict{String, DataFrame}()
    if ispath(path_sub_graph)
        file = jldopen(path_sub_graph, "r")
        for key in keys(file)
            subgraphs[key] = file[key]
        end
        close(file)
    else
        for node in male_nodes
            df_neighbors = df_male[(df_male[!, :"From Node ID"] .== node) .| (df_male[!, :"To Node ID"] .== node), :]
            subgraphs[node] = df_neighbors
        end
    end
    """get the dictionary of male, female edges and male and female matching
       male_edges (Dict): ("From Node ID", "To Node ") -> Weight
       female_edges (Dict): ("From Node ID", "To Node ") -> Weight
       matching (Dict): "Male Node ID" -> "Female Node ID"
       matching_female (Dict): "Female Node ID" -> "Male Node ID"
    """
    male_edges = Dict(zip(zip(df_male."From Node ID", df_male."To Node ID"), df_male.Weight))
    female_edges = Dict(zip(zip(df_female."From Node ID", df_female."To Node ID"), df_female.Weight))
    matching = Dict(zip(df_sol."Male Node ID", df_sol."Female Node ID"))
    matching_female = Dict(zip(df_sol."Female Node ID", df_sol."Male Node ID"))
    # initialize current score, new matching, current time, total swap steps
    cur_score = calculate_final_score(male_edges, female_edges, matching)
    new_matching = copy(matching)
    cur_time = now()
    steps = 0
    # initialize WandB config logger
    config = Dict(
        "algo" => "greedy_julia (random male node)", "init_score" => cur_score,
        "seed" => args["seed"]) 
    Wandb.login()
    lg = WandbLogger(project = "vnc_matching",
                 name = now_time,
                 mode="offline",
                 config = config)
    global_logger(lg)
    println("===============The random greedy search process is running=============")
    for iter = 1:args["max_iter"]
        # for every iter, select a random target male node and get its sub graph and female mapped graph 
        target_node = rand(male_nodes)
        df_target_subset = subgraphs[target_node]
        df_target_subset_mapped = map_female(df_target_subset, female_edges, new_matching)
        # initialize the best score delta and best matching
        best_delta = 0
        best_matching = Dict{String, String}()
        # for every selected male node, loop all the female nodes and swap its current matched female node with all other female nodes
        for female_node in female_nodes
            steps += 1
            # get the partner male node of current female node
            partner_node = matching_female[female_node]
            # calculate score delta and get the new matching
            delta_temp, new_matching_temp = cal_delta_random_greedy(
                partner_node, target_node, df_target_subset,
                df_target_subset_mapped, female_edges, subgraphs,
                new_matching)
            # if better delta, update best delta and best matching
            if delta_temp > best_delta
                best_delta = delta_temp
                best_matching = new_matching_temp
            end
        end
        # log at each male node
        time_elapsed = now() - cur_time
        Wandb.log(lg, Dict("score_increase" => best_delta, 
                            "cur_score" => Int(cur_score+best_delta),
                            "n_male_nodes" => iter, 
                            "time_elapsed" => time_elapsed, 
                            "steps" => steps))
        if best_delta > 0
            # update score and matching, save new matching to local
            filename_rm = joinpath("..", "sol", "$(Int(cur_score))_$(now_time).csv")
            if ispath(filename_rm)
                rm(filename_rm)
            end
            cur_score += best_delta
            new_matching = best_matching
            df_matching = transfer_df_matching(new_matching)
            save_sol(df_matching, cur_score, now_time)
            cur_time = now()
        end
    end
end
main()