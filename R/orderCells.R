# assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, last_pseudotime, curr_cell)
# {
#     nei <- NULL
#
#     cell_tree <- ordering_tree_res$subtree
#     curr_cell_pseudotime <- last_pseudotime
#     V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
#     V(cell_tree)[curr_cell]$parent =  V(cell_tree)[ nei(curr_cell, mode="in") ]$name
#     #print (curr_cell_pseudotime)
#
#     ordering_tree_res$subtree <- cell_tree
#     children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
#
#     for (child in children)	{
#         next_node <- V(cell_tree)[child]$name
#         delta_pseudotime <- dist_matrix[curr_cell, next_node]
#         ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, dist_matrix, last_pseudotime + delta_pseudotime, next_node)
#     }
#
#     return (ordering_tree_res)
# }
#res <- assign_pseudotime_helper(res, dist_matrix, 0.0, res$root)

curr_state <- 1

stree <- DDRTree_res$stree +  t(DDRTree_res$stree) != 0
stree <- as.matrix(stree)
stree[stree = T] <- 1
stree[stree = 0] <- 0

#stree[upper.tri(stree)] <- 0
dimnames(stree) <- list(as.character(1:487), as.character(1:487))
stree_g <- graph.adjacency(stree, mode = "directed", diag = F, weighted = F)

res <- list(subtree = stree_g, root = "358")

assign_cell_state_helper <- function(ordering_tree_res, curr_cell)
{
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state = curr_state

    children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
    ordering_tree_res$subtree <- cell_tree
    message('curr_cell: ', curr_cell)
    message('children: ', children)

    if (length(children) == 1){
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name)
    }else{
        for (child in children)	{

            curr_state <<- curr_state + 1
            ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name)
        }
    }
    return (ordering_tree_res)
}

res <- assign_cell_state_helper(res, res$root)



