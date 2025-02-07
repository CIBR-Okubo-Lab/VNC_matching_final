include("./MatchingFunctions.jl")
using ArgParse
using CSV
using DataFrames
using DataFramesMeta
using Dates
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

        "--tol"
        help = "Bound for convergence, when convergence change the random swap direction to male/female"
        arg_type = Float16
        default = 5e-6
    end
    # init the random seed and file path
    args = parse_args(parser)
    Random.seed!(args["seed"])
    path_male = joinpath("..", "data", "male_connectome_graph.csv")
    path_female = joinpath("..", "data", "female_connectome_graph.csv")
    path_sol = joinpath("..", "sol", args["init_sol"])
    # read the male, female connection data and initial solution file
    df_male, df_female, df_sol = read_data(path_male, path_female, path_sol)
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

    # get the list of unique male and female edges
    unique_edges_female = collect(keys(female_edges))
    unique_edges_male = collect(keys(male_edges))
    # initialize the current score, Convergence tolerance and matching dictionary
    cur_score = calculate_final_score(male_edges, female_edges, matching)
    convergence = false    # 0(false) -> male, 1(true) -> female
    tol = args["tol"]
    new_matching = copy(matching)
    new_matching_female = copy(matching_female)

    # config the Wandb Logger
    Wandb.login()
    config = Dict("algo" => "Random Swap M/F Auto Version",
                "init_score" => cur_score,
                "seed" => args["seed"])
    lg = WandbLogger(project = "vnc_matching",
                 name = now_time,
                 mode="offline",
                 config = config)
    global_logger(lg)
    steps = 0   # count the total number of swap
    sum_delta = 0   # record the total sum of score_increase
    # main loop
    while true     
        steps += 1
        if convergence      # for female round
            # randomly choose a female edge
            random_female_edge = rand(unique_edges_female)
            # get the corresponding male edge according to current matching
            male_edge = (new_matching_female[random_female_edge[1]], new_matching_female[random_female_edge[2]])
            # swap the two male nodes in the corresponding male edge and calculate score delta 
            delta, new_matching_temp = cal_delta_random_swap(male_edge, df_male, female_edges, new_matching)
            # if score delta > 0, update current information and save new solution
            if delta > 0
                # delete the current solution file before saving new one
                filename_rm = joinpath("sol", "$(Int(cur_score))_$(now_time).csv")
                if ispath(filename_rm)
                    rm(filename_rm)
                end
                # update current score and male and female matching
                cur_score = cur_score + delta
                new_matching = new_matching_temp 
                new_matching_female[random_female_edge[1]], new_matching_female[random_female_edge[2]] = new_matching_female[random_female_edge[2]], new_matching_female[random_female_edge[1]]
                # converting dictionary type matching to DataFrame and save it
                df_matching = transfer_df_matching(new_matching)
                save_sol(df_matching, cur_score, now_time)
                sum_delta += delta
            # if score delta <= 0, do nothing
            else
                delta = 0
            end
        else        # for male round
            # randomly choose a male edge(pair of male nodes)
            random_male_edge = rand(unique_edges_male) 
            # swap the two male nodes in the selected male edge and calculate score delta 
            delta, new_matching_temp = cal_delta_random_swap(random_male_edge, df_male, female_edges, new_matching)
            # if score delta > 0, update current score and matching and save new solution
            if delta > 0 
                # delete current solution file before saving new one 
                filename_rm = joinpath("sol", "$(Int(cur_score))_$(now_time).csv")
                if ispath(filename_rm)
                    rm(filename_rm)
                end
                cur_score = cur_score + delta
                new_matching = new_matching_temp # updated matching
                df_matching = transfer_df_matching(new_matching)
                save_sol(df_matching, cur_score, now_time)
                sum_delta += delta
            else    # no score increase
                delta = 0
            end
        end
        # log the score increase, current score and total swapping steps in Wandb
        Wandb.log(lg, 
            Dict("score_increase" => delta, 
                "cur_score" => Int(cur_score),
                "steps" => steps))
        # check convergence every One million steps
        if steps % 1e6 == 0
            if sum_delta / 1e6 < tol
                convergence = ~convergence
                println("change of direction of swap")
            end
            sum_delta = 0
        end
    end
    close(lg)
end
main()