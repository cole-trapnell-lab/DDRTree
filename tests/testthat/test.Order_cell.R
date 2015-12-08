#compare the results with the original ordering and pseudotime assignment
reduceDimension <- function (cds, max_components = 2, use_irlba = TRUE, pseudo_expr = 1,
          batch = NULL, batch2 = NULL, covariates = NULL, use_vst = FALSE,
          verbose = FALSE, ...)
{
    FM <- exprs(cds)
    if (is.null(use_vst) && cds@expressionFamily@vfamily == "negbinomial") {
        use_vst = TRUE
        pseudo_expr = 0
    }
    if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial") {
        checkSizeFactors(cds)
        size_factors <- sizeFactors(cds)
        FM <- t(t(FM)/size_factors)
    }
    if (is.null(fData(cds)$use_for_ordering) == FALSE && nrow(subset(fData(cds),
                                                                     use_for_ordering == TRUE)) > 0)
        FM <- FM[fData(cds)$use_for_ordering, ]
    if (cds@expressionFamily@vfamily == "binomialff") {
        ncounts <- FM
        ncounts[ncounts != 0] <- 1
        FM <- t(t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
    }
    if (cds@expressionFamily@vfamily != "binomialff") {
        FM <- FM + pseudo_expr
    }
    FM <- FM[matrixStats::rowSds(FM) > 0, ]
    if (cds@expressionFamily@vfamily != "binomialff") {
        if (use_vst) {
            VST_FM <- vstExprs(cds, round_vals = FALSE)
            if (is.null(VST_FM) == FALSE) {
                FM <- VST_FM
                FM <- FM[fData(cds)$use_for_ordering, ]
            }
            else {
                stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
            }
        }
        else {
            FM <- log2(FM)
        }
    }
    if (is.null(batch) == FALSE || is.null(batch2) == FALSE ||
        is.null(covariates) == FALSE) {
        if (verbose)
            message("Removing batch effects")
        FM <- limma::removeBatchEffect(FM, batch = batch, batch2 = batch2,
                                       covariates = covariates)
        if (cds@expressionFamily@vfamily != "binomialff") {
            if (use_vst == FALSE) {
                FM <- 2^FM
            }
        }
    }
    if (verbose)
        message("Reducing to independent components")
    init_ICA <- ica_helper(t(FM), max_components, use_irlba = use_irlba,
                           ...)
    x_pca <- t(t(FM) %*% init_ICA$K)
    W <- t(init_ICA$W)
    weights <- W
    A <- t(solve(weights) %*% t(init_ICA$K))
    colnames(A) <- colnames(weights)
    rownames(A) <- rownames(FM)
    S <- weights %*% x_pca
    rownames(S) <- colnames(weights)
    colnames(S) <- colnames(FM)
    reducedDimW(cds) <- W
    reducedDimA(cds) <- A
    reducedDimS(cds) <- S
    reducedDimK(cds) <- init_ICA$K
    cds
}

orderCells <- function (cds, num_paths = 1, reverse = FALSE, root_cell = NULL,
          scale_pseudotime = F)
{
    adjusted_S <- t(cds@reducedDimS)
    dp <- as.matrix(dist(adjusted_S))
    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cc_ordering <- extract_good_branched_ordering(res$subtree,
                                                  res$root, cellPairwiseDistances(cds), num_paths, reverse)
    row.names(cc_ordering) <- cc_ordering$sample_name
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),
                                         ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), ]$cell_state
    pData(cds)$Parent <- cc_ordering[row.names(pData(cds)), ]$parent
    if (scale_pseudotime) {
        cds <- scale_pseudotime(cds)
    }
    cds
}

plot_spanning_tree <- function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE,
          show_backbone = TRUE, backbone_color = "black", markers = NULL,
          show_cell_names = FALSE, cell_size = 1.5, cell_link_size = 0.75,
          cell_name_size = 2, show_all_lineages = F)
{
    gene_short_name <- NULL
    sample_name <- NULL
    lib_info_with_pseudo <- pData(cds)
    S_matrix <- reducedDimS(cds)
    if (is.null(S_matrix)) {
        stop("You must first call reduceDimension() before using this function")
    }
    ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
    colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
    ica_space_df$sample_name <- row.names(ica_space_df)
    ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo,
                                     by.x = "sample_name", by.y = "row.names")
    dp_mst <- minSpanningTree(cds)
    if (is.null(dp_mst)) {
        stop("You must first call orderCells() before using this function")
    }
    edge_list <- as.data.frame(get.edgelist(dp_mst))
    colnames(edge_list) <- c("source", "target")
    edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name",
                     by.y = "source", all = TRUE)
    edge_df <- plyr::rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1",
                                       ICA_dim_2 = "source_ICA_dim_2"))
    edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name",
                                                          "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name",
                     all = TRUE)
    edge_df <- plyr::rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1",
                                       ICA_dim_2 = "target_ICA_dim_2"))
    diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst,
                                                           weights = NA)]$name))
    colnames(diam) <- c("sample_name")
    diam <- plyr::arrange(merge(ica_space_with_state_df, diam,
                                by.x = "sample_name", by.y = "sample_name"), Pseudotime)
    if (show_all_lineages) {
        pro_state_pseudotime <- diam[as.numeric(diam$State) ==
                                         min(as.numeric(diam$State)), "Pseudotime"]
        bifurcation_sample <- diam[which(diam$Pseudotime == max(pro_state_pseudotime)),
                                   ]
        bifurcation_sample$State <- diam$State[which(diam$Pseudotime ==
                                                         max(pro_state_pseudotime)) + 1]
        diam <- rbind(diam[1:which(diam$Pseudotime == max(pro_state_pseudotime)),
                           ], bifurcation_sample, diam[(which(diam$Pseudotime ==
                                                                  max(pro_state_pseudotime)) + 1):nrow(diam), ])
        no_diam_states <- setdiff(lib_info_with_pseudo$State,
                                  lib_info_with_pseudo[diam[, 1], "State"])
        for (state in no_diam_states) {
            state_sample <- ica_space_with_state_df[ica_space_with_state_df$State ==
                                                        state, "sample_name"]
            subset_dp_mst <- induced.subgraph(dp_mst, state_sample,
                                              impl = "auto")
            subset_diam <- as.data.frame(as.vector(V(subset_dp_mst)[get.diameter(subset_dp_mst,
                                                                                 weights = NA)]$name))
            colnames(subset_diam) <- c("sample_name")
            subset_diam <- plyr::arrange(merge(ica_space_with_state_df,
                                               subset_diam, by.x = "sample_name", by.y = "sample_name"),
                                         Pseudotime)
            subset_diam$State <- state
            bifurcation_sample$State <- state
            diam <- rbind(diam, bifurcation_sample, subset_diam)
        }
    }
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in%
                                    markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(exprs(cds[row.names(markers_fData),
                                                      ]))
            markers_exprs <- merge(markers_exprs, markers_fData,
                                   by.x = "Var1", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) >
        0) {
        edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name",
                         by.y = "Var2")
        g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1,
                                        y = source_ICA_dim_2, size = log10(value + 0.1))) +
            facet_wrap(~feature_label)
    }
    else {
        g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1,
                                        y = source_ICA_dim_2))
    }
    if (show_tree) {
        g <- g + geom_segment(aes_string(xend = "target_ICA_dim_1",
                                         yend = "target_ICA_dim_2", color = color_by), size = 0.3,
                              linetype = "solid", na.rm = TRUE)
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) >
        0) {
        for (i in unique(diam[, "State"])) {
            g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2),
                               color = I(backbone_color), size = I(cell_link_size),
                               data = subset(diam, State == i), na.rm = TRUE)
        }
    }
    else {
        g <- g + geom_point(aes_string(color = color_by), size = I(cell_size),
                            na.rm = TRUE)
    }
    if (show_backbone) {
        if (backbone_color %in% colnames(diam))
            g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2),
                               color = diam[, backbone_color], size = I(cell_link_size),
                               data = diam, na.rm = TRUE)
        else g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2),
                                color = I(backbone_color), size = I(cell_link_size),
                                data = diam, na.rm = TRUE)
    }
    if (show_cell_names) {
        g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
    }
    g <- g + theme(panel.border = element_blank(), axis.line = element_line()) +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
        ylab("Component 1") + xlab("Component 2") + theme(legend.position = "top",
                                                          legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill = "white"))
    g
}

