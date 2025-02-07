module MatchingFunctions
using CSV
using DataFrames
using DataFramesMeta
using Statistics
# export the functions 
export read_data, calculate_final_score, map_female, cal_score, transfer_df_matching, save_sol, cal_delta_random_greedy, cal_delta_random_swap

function read_data(path_male, path_female, path_sol)
    """
    Reads data from three CSV files containing male and female connectome graphs, and benchmark matching data

    Args:
        path_male (String): path to the CSV file containing male connectome graph data
        path_female (String): path to the CSV file containing female connectome graph data
        path_benchmark (String): path to the CSV file containing benchmark matching data

    Returns:
        tuple: a tuple containing three pandas DataFrames: df_male, df_female, and df_sol
    """
    df_male = CSV.read(path_male, DataFrame)
    df_female = CSV.read(path_female, DataFrame)
    df_sol = CSV.read(path_sol, DataFrame)
    return df_male, df_female, df_sol
end;

function calculate_final_score(male_edges::Dict, female_edges::Dict, matching::Dict)
    """
    Organizer's code to calculate alignment score

    Parameters:
        male_edges (Dict): male edges information. keys: tuples of male edges, values: edge weights
        female_edges (Dict): female edges information. keys: tuples of female edges, values: edge weights
        matching (Dict): male to female matching. keys: male nodes, values: female nodes

    Returns:
        alignment (Int): alignment score
    """
    alignment = 0
    for (_, item) in enumerate(male_edges)
        edge_weight = item[2]
        female_nodes = (matching[item[1][1]], matching[item[1][2]])
        alignment += min(edge_weight, get(female_edges, female_nodes, 0))
    end
    return alignment
end;

function map_female(df::DataFrame, female_edges::Dict, matching::Dict)
    """
    Map the corresponding female nodes to the given dataframe

    Parameters:
        df (DataFrame): dataframe to map the female nodes
        female_edges (Dict): female edges information. keys: tuples of female edges, values: edge weights
        matching (Dict): male to female matching. keys: male nodes, values: female nodes

    Returns:
        df (DataFrame): the dataframe with the female nodes mapped in
                      'From Node ID Female', 'To Node ID Female', and 'Weight Female'
    """
    df = copy(df)
    df[!, :"From Node ID Female"] = [matching[id] for id in df."From Node ID"]
    df[!, :"To Node ID Female"] = [matching[id] for id in df."To Node ID"]
    df[!, :"Weight Female"] = [Int(get(female_edges, (row[:"From Node ID Female"], row[:"To Node ID Female"]), 0.0)) for row in eachrow(df)]
    return df
end;

function cal_score(df::DataFrame)
    """
    Calculate the alignment score given a dataframe of a subset of male nodes' neighbors

    Parameters:
        df (DataFrame): dataframe of a subset of male nodes' neighbors
            containing columns 'From Node ID', 'To Node ID', 'Weight',
            'From Node ID Female', 'To Node ID Female', 'Weight Female'

    Returns:
        score (Int): alignment score
    """
    score = sum(min.(df[!, :"Weight Female"], df[!, :"Weight"]))
    return score
end;

function transfer_df_matching(matching::Dict)
    """
    Converting dictionary type matching to DataFrame for better saving

    Parameters:
        matching (Dict): male to female matching. keys: male nodes, values: female nodes

    Returns:
        df (DataFrame): two columns, "Male Node ID", "Female Node ID"
    """
    df = DataFrame(
        "Male Node ID" => String[], 
        "Female Node ID" => String[]
    )
    for (i,j) in zip(keys(matching), values(matching))
        push!(df, (i,j))
    end
    return df
end;

function save_sol(df_matching::DataFrame, cur_score::Int, now_time::String)
    """
    Save the solution to a CSV file.

    Parameters:
        matching (Dict): male to female matching. keys: male nodes, values: female nodes
        cur_score (Int): alignment score
        now_time (String): current time
    """
    filename = joinpath("sol", "$(Int(cur_score))_$(now_time).csv")
    CSV.write(filename, df_matching)
end;

function cal_delta_random_swap(nodes_pair::Tuple, df_male::DataFrame, female_edges::Dict, matching::Dict)
    """
    Calculate the increase in alignment score after swapping two male nodes mapping

    Parameters:
        nodes_pair (Tuple): two male nodes to swap
        df_male (DataFrame): dataframe of male connectome
        female_edges (Dict): female edges information. keys: tuples of female edges, values: edge weights
        matching (Dict): male to female matching. keys: male nodes, values: female nodes

    Returns:
        tuple: (delta, new_matching)
            delta (Int): increase in alignment score
            new_matching (Dict): the new male to female matching after swapping two male nodes mapping
    """
    # create a subset of 2 male nodes' neighbors
    df_neighbors = @subset(df_male, in.(:"From Node ID", Ref(nodes_pair)) .|| in.(:"To Node ID", Ref(nodes_pair)))

    # map the corresponding female nodes
    df_matching_prev = map_female(df_neighbors, female_edges, matching)

    # calculate previous matching score
    score_prev = cal_score(df_matching_prev)

    # create a new matching after swapping the nodes
    new_matching = copy(matching)
    new_matching[nodes_pair[1]], new_matching[nodes_pair[2]] = new_matching[nodes_pair[2]], new_matching[nodes_pair[1]]

    # map the corresponding female nodes after swapping
    df_matching = map_female(df_neighbors, female_edges, new_matching)

    # calculate the new matching score
    score = cal_score(df_matching)
    delta = score - score_prev
    return delta, new_matching
end;

function cal_delta_random_greedy(partner_node, target_node, df_target_subset::DataFrame,
    df_target_subset_mapped::DataFrame, female_edges::Dict, subgraphs::Dict,
    matching::Dict)
    """
    Calculate the change in alignment score after swapping the nodes in the matching

    Parameters:
        partner_node (String): the partner node to be swapped
        target_node (String): the target node to be swapped
        df_target_subset (DataFrame): dataframe subset of the target node's neighbors 
        df_target_subset_mapped (DataFrame): dataframe with the target node's neighbors mapped to female nodes 
        female_edges (Dict): female edges with tuples as keys and edge weights as values 
        subgraphs (Dict): subgraphs with node IDs as keys and DataFrames as values 
        matching (Dict): the current male to female node matching 

    Returns:
        tuple: a tuple containing:
            - delta (Int): the change in alignment score after the swap 
            - temp_matching (Dict): the updated matching dictionary after the node swap 
    """
    # get the sub_graph of partner node and map the sub_graph with the female edges
    df_partner_subset = subgraphs[partner_node]
    df_partner_subset_mapped = map_female(df_partner_subset, female_edges, matching)
    # calculate the current sub_graph score
    df_mapped = vcat(df_target_subset_mapped, df_partner_subset_mapped) |> unique
    prev_score = cal_score(df_mapped)
    # swap the partner node and target node in matching
    temp_matching = copy(matching)
    temp_matching[target_node], temp_matching[partner_node] = temp_matching[partner_node], temp_matching[target_node]
    # get the mapped results according to new temp_matching
    df_target_subset_mapped_temp = map_female(df_target_subset, female_edges, temp_matching)
    df_partner_subset_mapped_temp = map_female(df_partner_subset, female_edges, temp_matching)
    # calculate the sub_graph score after swapping
    df_mapped_temp = vcat(df_target_subset_mapped_temp, df_partner_subset_mapped_temp) |> unique
    cur_score = cal_score(df_mapped_temp)
    # calculate score delta after and before swapping
    delta = cur_score - prev_score
    return delta, temp_matching
end;

end