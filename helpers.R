# helpers.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 06/09/2023
# updated 09/13/2023
# Description: This script creates several functions used to analyze CRISPR gene effect data from the DepMap dataset
#
# For this first function, what we're interested in are mutational features that are enriched in the cell lines most dependent on the genes of interest. To find these, we'll first rank order cell lines based on the CRISPR screening knock-out effect then perform mutation enrichment analysis across the gene effect distribution by using "gene set enrichment" analysis. We will then perform a similar analysis at the cancer type level in order to identify cancer type enrichment based on the gene effect distribution.

crispr_dependency_enrichment <-
  function(genes_of_interest,
           cancer_type) {
    # subset to only genes of interest and information on the type of cancer the cell lines represent
    sub_data = effect_scores |>
      select(
        ModelID,
        c(genes_of_interest),
        StrippedCellLineName,
        DepmapModelType,
        OncotreeSubtype,
        OncotreePrimaryDisease,
        OncotreeLineage
      ) |>
      unique()
    
    # subset the data to the desired cancer of interest
    if (cancer_type %ni% c("PAN-CANCER", "PAN CANCER", "Pan-Cancer", "Pan Cancer", "Pan-cancer", "pan-cancer", "pan cancer")) {
      sub_data <- sub_data |>
        subset(DepmapModelType == cancer_type)
    }
    
    # in the current iteration of this function, only one input gene will be used, however this for loop is intended to work for lists of genes in the future
    for (i in 1:length(genes_of_interest)) {
      # create lists of plots and data in order to plot and download on the shiny app interface
      plot_list <- list()
      results_list <- list()
      
      # find index for gene of interest
      gene_column <-
        which(colnames(sub_data) == genes_of_interest[i])
      
      # rank-order the cell lines based on the gene effect score for the gene of interest
      sub_data_ranked <- sub_data |>
        arrange_at(which(colnames(sub_data) == genes_of_interest[i]))  |>
        unique() |>
        mutate(gene_dependency_order = row_number(),)
      
      # rename column names for easier data manipulation
      colnames(sub_data_ranked) <-
        c(
          "ModelID",
          "gene",
          "StrippedCellLineName",
          "DepmapModelType",
          "OncotreeSubtype",
          "OncotreePrimaryDisease",
          "OncotreeLineage",
          "gene_dependency_order"
        )
      
      sub_data_ranked <- sub_data_ranked |>
        mutate(dependency_group = case_when(gene < 0 ~ "More_dependent",
                                            TRUE ~ "Less_dependent"))
      
      # plot the dependency distribution for the given gene of interest 
      # gene effect distribution ----
      p1 <-
        ggplot(sub_data_ranked, aes(reorder(ModelID,gene), gene)) +
        geom_bar(aes(fill = gene), stat = "identity") +
        scale_color_viridis(
          option = "viridis",
          direction = -1,
          name = paste(genes_of_interest[i], "effect score", sep = "")
        ) +
        scale_fill_viridis(
          option = "viridis",
          direction = 1,
          name = paste(genes_of_interest[i], "effect score", sep = " ")
        ) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
        ) +
        labs(
          title = paste(
            cancer_type,
            "\n",
            genes_of_interest[i],
            "effect distribution",
            sep = " "
          ),
          x = "Cell line rank order\n<- more dependent lines | less dependent lines ->",
          y = "Gene Effect Score"
        )
      
      # now indicate the tissue of origin for the different cell lines plotted
      # gene effect distribution by cancer ----
      p2 <-
        ggplot(sub_data_ranked, aes(reorder(ModelID, gene), gene)) +
        geom_bar(aes(fill = OncotreeLineage), stat = "identity") +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title.align = 0.5
        ) +
        labs(
          title = paste(
            cancer_type,
            "\n",
            genes_of_interest[i],
            "effect distribution",
            sep = " "
          ),
          x = "Cell line rank order\n<- less dependent lines | more dependent lines ->",
          y = "Gene Effect Score"
        ) +
        guides(fill = guide_legend(title = "Tissue"))
      
      # save data to list
      results_list[[1]] <- sub_data_ranked
      
      # add plots to plot list
      plot_list[[1]] <- p1
      plot_list[[2]] <- p2
      
      # Mutation-specific enrichment ----
      # add the mutation data to the gene effect ranked data
      sub_data_ranked_muts <-
        left_join(sub_data_ranked, mutations, by = "ModelID")
      
      # select mutations that are present in multiple cell lines
      # need to specify the number of cell lines required to perform the analysis based on a pan-cancer or cancer-specific filter. these numbers are arbitrary but if desired, a percentage filter could be implimented here in the future
      # define the number based on the type of analysis performed
      n_cells <- if_else(cancer_type %in% c("PAN-CANCER", "PAN CANCER", "Pan-Cancer", "Pan Cancer", "Pan-cancer", "pan-cancer", "pan cancer"), 15, 5)
      
      # subset the data
      sub_data_ranked_muts <- sub_data_ranked_muts |>
        group_by(HugoSymbol) |>
        mutate(n_lines = n_distinct(ModelID)) |>
        ungroup() |>
        subset(n_lines >= n_cells)
      
      # Methodology to determine mutation enrichment across cell lines
      # To do this, we will perform enrichment analysis using fgsea and instead of sets of genes we will use sets of cell lines.
      # The pathway is the set of cell lines with a given mutation and the fold-change is the CRISPR score.
      # You then test whether that "pathway" of cell  lines is enriched in cells with high CRISPR score.
      # I think, in order to do this, what I need to do is, for all mutations making it into the enrichment part of the function,
      # generate the "pathways" file for all mutation sets
      # For this we will need rank-ordered list of samples and their group assignments, which will be based on mutation presence
      
      # rank-ordered list of samples and their group assignments
      ranked_samples <- sub_data_ranked_muts |>
        select(
          ModelID,
          which(colnames(sub_data) == genes_of_interest[i]),
          gene_dependency_order,
          dependency_group
        ) |>
        unique()
      
      # Mutations and their presence across the samples
      muts <-
        c(unique(sub_data_ranked_muts$HugoSymbol))
      
      # need to specify if the analysis is done using the pan-cancer cohort or cancer-specific subset
      # these values will be used in the gsea analyses performed downstream
      if (n_cells >= 15) {
        group_size = n_cells
        t_path <- 10
        b_path <- 10
        mut_label1 <- 11
        mut_label2 <- 11
      }
      if (n_cells < 15) {
        group_size = n_cells
        t_path <- 2
        b_path <- 2
        mut_label1 <- 3
        mut_label2 <- 2
      }
      
      # Function to determine mutation enrichment
      # this function takes in a rank-ordered list of cell lines in addition to a metadata list of group assignments
      mutation_enrichment <-
        function(ranked_samples, muts) {
          mutated_cell_line_pathways <- list()
          
          # Generate and add multiple lists to the main list
          for (i in 1:length(muts)) {
            # identify the mutated cell lines
            mutated_cell_line_pathways[muts[i]] <-
              sub_data_ranked_muts |>
              subset(HugoSymbol == muts[i]) |>
              select(ModelID) |>
              unique() |>
              `colnames<-`(muts[i]) |>
              c()
          }
          
          # randomize the ties
          genesTables_depmap <- ranked_samples |>
            na.omit() |>
            mutate(rank = rank(muts[i],  ties.method = "random")) |>
            dplyr::arrange(desc(rank))
          
          # change the column names in order to pass the correct column of values into the stats argument of the fgsea function
          colnames(genesTables_depmap) <-
            c("ModelID",
              "ko_gene",
              "gene_dependency_order",
              "dependency_group",
              "rank")
          
          depmap_Ranks = genesTables_depmap$ko_gene
          names(depmap_Ranks) = c(genesTables_depmap$ModelID)
          
          # run fgsea
          fgseaResult_depmap <-
            fgsea(
              pathways = mutated_cell_line_pathways,
              stats = depmap_Ranks,
              minSize = group_size,
              maxSize = 500
            )
          
          # select the "pathways" or in this case genes that show the most differential enrichment
          topPathwaysUp <-
            fgseaResult_depmap[ES > 0][head(order(pval), n = t_path), pathway]
          
          topPathwaysDown <-
            fgseaResult_depmap[order(ES)][head(order(ES), n = b_path), pathway]
          
          topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
          
          collapsedPathways <-
            collapsePathways(fgseaResult_depmap[order(pval)][padj < 0.01],
                             mutated_cell_line_pathways,
                             depmap_Ranks)
          mainPathways_depmap <-
            fgseaResult_depmap[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
          
          fgseaResult_depmap_final <-
            as.data.frame(fgseaResult_depmap)
          
          # annotate the top and bottom genes with the most differentially enriched scores
          fgseaResult_depmap_final <-
            fgseaResult_depmap_final |>
            arrange_at("NES") |>
            mutate(
              rank = row_number(),
              mutation_label = case_when(
                rank < mut_label1 |
                  rank > (nrow(fgseaResult_depmap_final) - mut_label2) ~ pathway,
                TRUE ~ NA
              )
            )
          return(
            list(
              fgseaResult_depmap_final = data.table(fgseaResult_depmap_final),
              mutated_cell_line_pathways = mutated_cell_line_pathways,
              depmap_Ranks = depmap_Ranks,
              topPathways = topPathways
            )
          )
        }
      
      # Call the function to determine mutation enrichment
      mutation_enrichment <-
        mutation_enrichment(ranked_samples, muts)
      
      # plot the differential enrichment scores for the given gene of interest
      # differential enrichment scores ----
      p3 <-
        ggplot(
          mutation_enrichment$fgseaResult_depmap_final,
          aes(
            color = NES,
            reorder(pathway, -NES),
            NES,
            label = mutation_label
          )
        ) +
        geom_point() +
        geom_bar(aes(fill = NES), stat = "identity") +
        scale_color_viridis(name = "Normalized\nenrichment\nscore", direction = -1) +
        scale_fill_viridis(direction = -1) +
        geom_label_repel(inherit.aes = TRUE, max.overlaps = 15) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title.align = 0.5
        ) +
        labs(
          title = paste(
            "Distribution of mutation enrichment\nbased on",
            genes_of_interest[i],
            "dependency"
          ),
          x = "Mutation rank order\n<- enriched in more dependent lines | enriched in less dependent lines ->",
          y = "Mutation Normalized Enrichment Score"
        ) +
        guides(fill = FALSE)
      
      # save plot and raw data
      results_list[[2]] <- as.data.frame(mutation_enrichment$fgseaResult_depmap_final)
      
      # add plot to plot list
      plot_list[[3]] <- p3
      
      # plot most differentially enriched mutations on their own
      plot_data_depmap = as.data.frame(mutation_enrichment$fgseaResult_depmap_final) |>
        subset(pathway %in% mutation_enrichment$topPathways)
      
      plot_data_depmap$mutation_label <- plot_data_depmap$pathway
      
      plot_data_depmap <- na.omit(plot_data_depmap)
      
      plot_data_depmap = plot_data_depmap |>
        mutate(
          `Enriched in` = case_when(
            NES < 0 & row_number() < 11 ~ "Less dependent lines",
            NES > 0 &
              row_number() > 10 ~ "More dependent lines"
          ),
          significance = case_when(padj < 0.1 ~ "padj < 0.1",
                                   padj > 0.1 ~ "padj > 0.1")
        )
      
      plot_data_depmap$significance = factor(plot_data_depmap$significance,
                                             levels = c("padj > 0.1", "padj < 0.1"))
      
      # dynamically set plotting aesthetics
      y_limit <- nrow(plot_data_depmap) + 1
      left_arrow <- min(plot_data_depmap$NES)
      right_arrow <- max(plot_data_depmap$NES)
      left_text <- (left_arrow - .1) / 2
      right_text <- (right_arrow + .1) / 2
      
      # plot top genes ----
      p4 <-
        ggplot(
          plot_data_depmap |>
            arrange(pathway, NES),
          aes(
            x = NES,
            y = reorder(pathway, -NES),
            fill = `Enriched in`,
            alpha = sort(significance, increasing = T)
          )
        ) +
        geom_point(color = "black", shape = 21, aes(size = size)) +
        scale_fill_manual(values = c(
          "Less dependent lines" = "#fde725",
          "More dependent lines" = "#440154"
        )) +
        theme_cowplot() +
        ylab(NULL) +
        coord_cartesian(ylim = c(0, y_limit), clip = "off") +
        geom_vline(
          xintercept = 0,
          color = "#969696",
          linetype = "dashed",
          size = .5
        ) +
        scale_alpha_discrete(range = c(0.5, 1), guide = FALSE) +
        annotate(
          "segment",
          x = 0.1,
          y = nrow(plot_data_depmap) + 1,
          xend = right_arrow,
          yend = nrow(plot_data_depmap) + 1,
          size = 5,
          linejoin = "mitre",
          arrow = arrow(type = "closed", length = unit(0.01, "npc")),
          color = "#440154"
        ) +
        annotate(
          "segment",
          x = -0.1,
          y = nrow(plot_data_depmap) + 1,
          xend = left_arrow,
          yend = nrow(plot_data_depmap) + 1,
          size = 5,
          linejoin = "mitre",
          arrow = arrow(type = "closed", length = unit(0.01, "npc")),
          color = "#fde725"
        ) +
        annotate(
          "text",
          x = right_text,
          y = nrow(plot_data_depmap) + 1,
          label = "More dependent",
          color = "white"
        ) +
        annotate(
          "text",
          x = left_text,
          y = nrow(plot_data_depmap) + 1,
          label = "Less dependent",
          color = "black"
        ) +
        guides(
          fill = guide_legend(
            override.aes = list(size = 3),
            title = "Enriched in:",
            order = 2
          ),
          size = guide_legend(order = 1, title = "Number of\ncell lines"),
          alpha = guide_legend(override.aes = list(size = 3), title = "Significance")
        ) +
        ggtitle(
          paste(
            "Top differentially enriched mutations based on ",
            genes_of_interest[i],
            " dependency",
            sep = ""
          )
        ) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title.align = 0.5
        )
      
      # add data to data list
      results_list[[3]] <- plot_data_depmap
      
      # add plot to plot list
      plot_list[[4]] <- p4
      
      # plot the distribution of the top and bottom mutations across the dependency distribution for the gene of interest
      # select mutations for plotting
      muts_for_plotting1 <-
        mutation_enrichment$fgseaResult_depmap_final |>
        select(mutation_label, NES) |>
        unique() |>
        na.omit() |>
        mutate(
          direction = case_when(
            NES > 0 ~ "Enriched in less dependent lines",
            TRUE ~ "Enriched in more dependent lines"
          )
        ) |>
        `colnames<-`(c("variable", "NES", "direction"))
      
      # select only the top 5 mutations per group for plotting
      muts_for_plotting <- muts_for_plotting1 %>%
        group_by(direction) %>%
        mutate(rank = ifelse(direction > 0, rank(desc(NES)), rank(NES))) %>%
        filter(rank <= 5) %>%
        ungroup() %>%
        select(-rank)
      
      # sub_data_ranked_muts contains all mutations and effect scores
      sub_data_ranked_final <- sub_data_ranked |>
        select(ModelID, gene) |>
        tibble::add_column(!!!set_names(as.list(rep(
          NA, length(muts_for_plotting$variable)
        )), nm = muts_for_plotting$variable)) |>
        as_tibble()
      
      # find unique mutations
      mutations_2 <- c(unique(muts_for_plotting$variable))
      
      # Iterate through each mutation
      for (k in seq_along(mutations_2)) {
        # find cases with mutations
        mut_samples <- sub_data_ranked_muts |>
          subset(HugoSymbol == mutations_2[k]) |>
          select(ModelID) |>
          unique()
        
        # Separate the ranked samples into two groups based on mutation presence
        sub_data_ranked_final <- sub_data_ranked_final |>
          mutate(
            !!mutations_2[k] := case_when(ModelID %in% mut_samples$ModelID ~ "Mut",
                                          TRUE ~ "WT")
          )
      }
      sub_data_ranked_final_m <- sub_data_ranked_final |>
        select(-ModelID) |>
        reshape2::melt(id = "gene") |>
        mutate(rank = row_number())
      
      sub_data_ranked_final_m <-
        left_join(sub_data_ranked_final_m, muts_for_plotting, by = "variable")
      
      colnames(sub_data_ranked_final_m) <-
        c("gene",
          "variable",
          "Cell line mutation status",
          "rank",
          "NES",
          "direction")
      
      # plot mutated cell lines in gene effect distribution ----
      p5 <-
        ggplot(sub_data_ranked_final_m,
               aes(reorder(rank,gene),
                   gene)) +
        geom_bar(aes(fill = `Cell line mutation status`, color = `Cell line mutation status`),
                 stat = "identity") +
        scale_color_manual(values = c("WT" = "lightgrey", "Mut" = "#cb181d")) +
        scale_fill_manual(values = c("WT" = "white", "Mut" = "#cb181d")) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title.align = 0.5,
          strip.text = element_text(size = 10)
        ) +
        labs(
          title = paste(genes_of_interest[i], "dependency distribution", sep = " "),
          x = "Cell line rank order\n<- more dependent lines | less dependent lines ->",
          y = "Gene Effect Score"
        ) +
        facet_wrap(direction ~ variable, nrow = 2)
      
      # add data to data list
      results_list[[4]] <- sub_data_ranked_final_m
      
      # add plot to plot list
      plot_list[[5]] <- p5
      
      # now, plot the normalized enrichment curves for the top differentially enriched mutations
      # select only the top 10 mutations per group for plotting
      muts_for_plotting <- muts_for_plotting1 %>%
        group_by(direction) %>%
        mutate(rank = ifelse(direction > 0, rank(desc(NES)), rank(NES))) %>%
        filter(rank <= 10) %>%
        ungroup() %>%
        select(-rank)
      
      diff_genes <- c(unique(muts_for_plotting$variable))
      
      nes_plots <- list()
      
      # Initialize an empty data frame to store the extracted data
      combined_data <- data.frame()
      
      for (j in seq_along(diff_genes)) {
        nes_plot_data <-
          plotEnrichment(
            mutation_enrichment$mutated_cell_line_pathways[[diff_genes[j]]],
            mutation_enrichment$depmap_Ranks
          ) + labs(title = paste(diff_genes[j]))
        
        nes_plot_data <-
          nes_plot_data$data # extract the data from the plotEnrichment results
        
        nes_plot_data$gene <-
          diff_genes[j] # annotate the gene of interest
        
        # now, add the NES score calculated earlier
        nes <- mutation_enrichment$fgseaResult_depmap_final |>
          subset(pathway == diff_genes[j]) |>
          select(NES) |>
          unique() |>
          as.numeric() |>
          round(digits = 2)
        
        nes_plot_data$NES <- nes
        
        # now, add the padj score calculated earlier
        padj <- mutation_enrichment$fgseaResult_depmap_final |>
          subset(pathway == diff_genes[j]) |>
          select(padj) |>
          unique() |>
          as.numeric() |>
          round(digits = 2)
        
        nes_plot_data$padj <- padj
        
        # Bind the extracted data to the combined data frame
        combined_data <- rbind(combined_data, nes_plot_data)
      }
      
      combined_data$rank <- as.numeric(combined_data$rank)
      combined_data$ES <- as.numeric(combined_data$ES)
      
      combined_data <- combined_data |>
        arrange(-NES) |>
        group_by(gene) |>
        mutate(NES = replace(NES, row_number() > 1, NA)) |>
        ungroup()
      
      # Convert facet variable to factor with desired order
      combined_data$gene <-
        factor(combined_data$gene, levels = unique(combined_data$gene))
      
      # NES curves ----
      p6 <- ggplot(combined_data,
                   aes(color = rank,
                       x = rank,
                       y = ES)) +
        geom_hline(yintercept = 0) +
        geom_point(size = 2) +
        geom_line() +
        # geom_bar(aes(fill = y), stat = "identity", width = 0.01) +
        scale_color_viridis(name = "Mutated cell line\n gene effect rank", direction = 1) +
        scale_fill_viridis(direction = -1) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title.align = 0.5,
          strip.text = element_text(size = 15)
        ) +
        xlab("Cell line rank order\n<- more dependent lines | less dependent lines ->") +
        ylab("Enrichment Score") +
        ggtitle(
          paste(
            "Most differentially enriched mutations by",
            genes_of_interest[i],
            "effect score"
          )
        ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill = FALSE) +
        facet_wrap( ~ gene)
      
      # add NES scores to each plot
      nes_data <- combined_data |>
        subset(!is.na(NES))
      
      # find x and y ranges in the data in order to plot the NES scores in a specific portion of the plot
      # x_pos <- max(combined_data$x) * 0.8
      x_pos <- n_distinct(sub_data$StrippedCellLineName) * 0.8
      y_pos <- max(combined_data$y) * 0.75
      
      p6 <- p6 +
        geom_text(
          color = "black",
          size = 3.5,
          data = nes_data,
          mapping = aes(
            x = x_pos,
            y = y_pos,
            label = paste("NES = ", NES, "\n", "padj = ", padj, sep = "")
          )
        )
      
      # save data
      results_list[[5]] <- combined_data
      
      # add plot to plot list
      plot_list[[6]] <- p6
      
      
      # now, we perform the same analyses we just performed at the individual mutation level but instead of grouping cell lines by mutation status we group by cancer type
      # Cancer-specific enrichment ----
      if (cancer_type %in% c("PAN-CANCER", "PAN CANCER", "Pan-Cancer", "Pan Cancer", "Pan-cancer", "pan-cancer", "pan cancer")) {
        # select cancer types that are represented by multiple cell lines and select a rank-ordered list of samples and their group assignments
        ranked_samples_cancer <- sub_data_ranked |>
          group_by(OncotreeSubtype) |>
          mutate(n_lines = n_distinct(ModelID)) |>
          ungroup() |>
          subset(n_lines >= 10) |>
          select(ModelID,
                 gene,
                 DepmapModelType,
                 gene_dependency_order,
                 dependency_group) |>
          unique()
        
        # Cell lines and their presence across the samples
        cancer_type <-
          c(unique(ranked_samples_cancer$DepmapModelType))
        
        # Function to determine cancer type enrichment across effect distribution
        cancer_enrichment <-
          function(ranked_samples_cancer,
                   cancer_type) {
            cancer_cell_line_pathways <- list()
            
            # Generate and add multiple lists to the main list
            for (i in seq_along(cancer_type)) {
              # identify the groups of cancer cell lines
              cancer_cell_line_pathways[cancer_type[i]] <-
                ranked_samples_cancer |>
                subset(DepmapModelType == cancer_type[i]) |>
                select(ModelID) |>
                unique() |>
                `colnames<-`(cancer_type[i]) |>
                c()
            }
            
            # randomize the ties
            cancersTables_depmap <- ranked_samples_cancer |>
              na.omit() |>
              mutate(rank = rank(gene,  ties.method = "random")) |>
              dplyr::arrange(desc(rank))
            
            # change the column names in order to pass the correct column of values into the stats argument of the fgsea function
            colnames(cancersTables_depmap) <-
              c("ModelID",
                "ko_gene",
                "gene_dependency_order",
                "dependency_group",
                "rank")
            
            depmap_cancer_Ranks = cancersTables_depmap$ko_gene
            names(depmap_cancer_Ranks) = c(cancersTables_depmap$ModelID)
            
            # run fgsea
            fgseaResult_cancers_depmap <-
              fgsea(
                pathways = cancer_cell_line_pathways,
                stats = depmap_cancer_Ranks,
                minSize = 10,
                maxSize = 500,
                nPermSimple = 100000
              )
            
            # select the "pathways" or in this case genes that show the most differential enrichment
            topPathwaysUp <-
              fgseaResult_cancers_depmap[ES > 0][head(order(pval), n = 10), pathway]
            
            topPathwaysDown <-
              fgseaResult_cancers_depmap[order(ES)][head(order(ES), n = 10), pathway]
            
            topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
            
            collapsedPathways <-
              collapsePathways(
                fgseaResult_cancers_depmap[order(pval)][padj < 0.01],
                cancer_cell_line_pathways,
                depmap_cancer_Ranks
              )
            mainPathways_depmap <-
              fgseaResult_cancers_depmap[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
            
            fgseaResult_cancers_depmap_final <-
              as.data.frame(fgseaResult_cancers_depmap)
            
            # annotate the top and bottom genes with the most differentially enriched scores
            fgseaResult_cancers_depmap_final <-
              fgseaResult_cancers_depmap_final |>
              arrange_at("NES") |>
              mutate(
                rank = row_number(),
                cancer_label = case_when(rank < 11 |
                                           rank > (
                                             nrow(fgseaResult_cancers_depmap_final) - 10
                                           ) ~ pathway,
                                         TRUE ~ NA)
              )
            return(
              list(
                fgseaResult_cancers_depmap_final = data.table(fgseaResult_cancers_depmap_final),
                cancer_cell_line_pathways = cancer_cell_line_pathways,
                depmap_Ranks = depmap_cancer_Ranks,
                topPathways = topPathways
              )
            )
          }
        
        # Call the function to determine mutation enrichment
        cancer_enrichment <-
          cancer_enrichment(ranked_samples_cancer, cancer_type)
        
        # plot the differential enrichment scores for the given gene of interest
        # differential enrichment scores ----
        p7 <-
          ggplot(
            cancer_enrichment$fgseaResult_cancers_depmap_final,
            aes(
              color = NES,
              reorder(pathway, -NES),
              NES,
              label = cancer_label
            )
          ) +
          geom_point() +
          geom_bar(aes(fill = NES), stat = "identity") +
          scale_color_viridis(name = "Normalized\nenrichment\nscore", direction = -1) +
          scale_fill_viridis(direction = -1) +
          geom_label_repel(inherit.aes = TRUE, max.overlaps = 15) +
          theme(
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = 'white'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title.align = 0.5
          ) +
          labs(
            title = paste(
              "Distribution of cancer enrichment\nbased on",
              genes_of_interest[i],
              "gene effect"
            ),
            x = "Cancer rank order\n<- enriched in more dependent lines | enriched in less dependent lines ->",
            y = "Mutation Normalized Enrichment Score"
          ) +
          guides(fill = FALSE)
        
        # add data to data list
        results_list[[6]] <- as.data.frame(cancer_enrichment$fgseaResult_cancers_depmap_final)
        
        # add plot to plot list
        plot_list[[7]] <- p7
        
        # plot most differentially enriched mutations on their own
        plot_data_depmap = as.data.frame(cancer_enrichment$fgseaResult_cancers_depmap_final) |>
          subset(pathway %in% cancer_enrichment$topPathways)
        
        plot_data_depmap$mutation_label <- plot_data_depmap$pathway
        
        plot_data_depmap <- na.omit(plot_data_depmap)
        
        plot_data_depmap = plot_data_depmap |>
          mutate(
            `Enriched in` = case_when(
              NES < 0 & row_number() < 11 ~ "Less dependent lines",
              NES > 0 &
                row_number() > 10 ~ "More dependent lines"
            ),
            significance = case_when(padj < 0.1 ~ "padj < 0.1",
                                     padj > 0.1 ~ "padj > 0.1")
          )
        
        plot_data_depmap$significance = factor(plot_data_depmap$significance,
                                               levels = c("padj > 0.1", "padj < 0.1"))
        
        # dynamically set plotting aesthetics
        y_limit <- nrow(plot_data_depmap) + 1
        left_arrow <- min(plot_data_depmap$NES)
        right_arrow <- max(plot_data_depmap$NES)
        left_text <- (left_arrow - .1) / 2
        right_text <- (right_arrow + .1) / 2
        
        # plot top cancers ----
        p8 <-
          ggplot(
            plot_data_depmap |>
              arrange(pathway, NES),
            aes(
              x = NES,
              y = reorder(pathway, -NES),
              fill = `Enriched in`,
              alpha = sort(significance, increasing = T)
            )
          ) +
          geom_point(color = "black",
                     shape = 21,
                     aes(size = size)) +
          scale_fill_manual(values = c(
            "Less dependent lines" = "#fde725",
            "More dependent lines" = "#440154"
          )) +
          theme_cowplot() +
          ylab(NULL) +
          coord_cartesian(ylim = c(0, y_limit), clip = "off") +
          geom_vline(
            xintercept = 0,
            color = "#969696",
            linetype = "dashed",
            size = .5
          ) +
          scale_alpha_discrete(range = c(0.5, 1), guide = FALSE) +
          annotate(
            "segment",
            x = 0.1,
            y = nrow(plot_data_depmap) + 1,
            xend = right_arrow,
            yend = nrow(plot_data_depmap) + 1,
            size = 5,
            linejoin = "mitre",
            arrow = arrow(type = "closed", length = unit(0.01, "npc")),
            color = "#440154"
          ) +
          annotate(
            "segment",
            x = -0.1,
            y = nrow(plot_data_depmap) + 1,
            xend = left_arrow,
            yend = nrow(plot_data_depmap) + 1,
            size = 5,
            linejoin = "mitre",
            arrow = arrow(type = "closed", length = unit(0.01, "npc")),
            color = "#fde725"
          ) +
          annotate(
            "text",
            x = right_text,
            y = nrow(plot_data_depmap) + 1,
            label = "More dependent",
            color = "white"
          ) +
          annotate(
            "text",
            x = left_text,
            y = nrow(plot_data_depmap) + 1,
            label = "Less dependent",
            color = "black"
          ) +
          guides(
            fill = guide_legend(
              override.aes = list(size = 3),
              title = "Enriched in:",
              order = 2
            ),
            size = guide_legend(order = 1, title = "Number of\ncell lines"),
            alpha = guide_legend(override.aes = list(size = 3), title = "Significance")
          ) +
          ggtitle(
            paste(
              "Top differentially enriched cancers based on ",
              genes_of_interest[i],
              " dependency",
              sep = ""
            )
          ) +
          theme(
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = 'white'),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title.align = 0.5
          )
        
        # add data to data list
        results_list[[7]] <- plot_data_depmap
        
        # add plot to plot list
        plot_list[[8]] <- p8
        
        # plot the distribution of the top and bottom cancers across the dependency distribution for the gene of interest
        # select cancers for plotting
        cancers_for_plotting1 <-
          cancer_enrichment$fgseaResult_cancers_depmap_final |>
          select(cancer_label, NES) |>
          unique() |>
          na.omit() |>
          mutate(
            direction = case_when(
              NES > 0 ~ "Enriched in less dependent lines",
              TRUE ~ "Enriched in more dependent lines"
            )
          ) |>
          `colnames<-`(c("variable", "NES", "direction"))
        
        cancers_for_plotting <- cancers_for_plotting1 %>%
          group_by(direction) %>%
          mutate(rank = ifelse(direction > 0, rank(desc(NES)), rank(NES))) %>%
          filter(rank <= 5) %>%
          ungroup() %>%
          select(-rank)
        
        # sub_data_ranked_muts contains all mutations and effect scores
        sub_data_ranked_final <- sub_data_ranked |>
          select(ModelID, DepmapModelType, gene) |>
          tibble::add_column(!!!set_names(as.list(rep(
            NA, length(cancers_for_plotting$variable)
          )), nm = cancers_for_plotting$variable)) |>
          as_tibble()
        
        # find unique cancers
        cancers_2 <- c(unique(cancers_for_plotting$variable))
        
        # Iterate through each cancer
        for (k in seq_along(cancers_2)) {
          # find cases for cancer of interest
          cancer_samples <- sub_data_ranked_final |>
            subset(DepmapModelType == cancers_2[k]) |>
            select(DepmapModelType) |>
            unique()
          
          # Separate the ranked samples into two groups based on mutation presence
          sub_data_ranked_final <- sub_data_ranked_final |>
            mutate(
              !!cancers_2[k] := case_when(
                DepmapModelType %in% cancer_samples$DepmapModelType ~ "Yes",
                TRUE ~ "No"
              )
            )
        }
        sub_data_ranked_final_m <- sub_data_ranked_final |>
          select(-DepmapModelType) |>
          reshape2::melt(id = "gene") |>
          mutate(rank = row_number())
        
        sub_data_ranked_final_m <-
          left_join(sub_data_ranked_final_m, cancers_for_plotting, by = "variable")
        
        colnames(sub_data_ranked_final_m) <-
          c("gene",
            "variable",
            "Cancer of interest",
            "rank",
            "NES",
            "direction")
        
        # plot the dependency distribution for the given gene of interest and visualize the cell lines used for downstream analyses
        # plot mutated cancer types in gene effect distribution ----
        p9 <-
          ggplot(sub_data_ranked_final_m,
                 aes(reorder(rank,gene),
                     gene)) +
          geom_bar(aes(fill = `Cancer of interest`, color = `Cancer of interest`),
                   stat = "identity") +
          scale_color_manual(values = c("No" = "lightgrey", "Yes" = "#cb181d")) +
          scale_fill_manual(values = c("No" = "lightgrey", "Yes" = "#cb181d")) +
          theme(
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = 'white'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title.align = 0.5,
            strip.text = element_text(size = 10)
          ) +
          labs(
            title = paste(genes_of_interest[i], "dependency distribution", sep = " "),
            x = "Cell line rank order\n<- more dependent lines | less dependent lines ->",
            y = "Gene Effect Score"
          ) +
          facet_wrap(direction ~ variable, nrow = 2)
        
        # add data to data list
        results_list[[8]] <- sub_data_ranked_final_m
        
        # add plot to plot list
        plot_list[[9]] <- p9
        
        # NES curves per cancer ----
        # now, plot the normalized enrichment curves for the top differentially enriched cancers
        diff_cancers <-
          unique(cancer_enrichment$fgseaResult_cancers_depmap_final$cancer_label) |>
          na.omit() |>
          c()
        
        nes_plots <- list()

        # Initialize an empty data frame to store the extracted data
        combined_data2 <- data.frame()
        
        for (j in seq_along(diff_cancers)) {
          nes_plot_data <-
            plotEnrichment(
              cancer_enrichment$cancer_cell_line_pathways[[diff_cancers[j]]],
              cancer_enrichment$depmap_Ranks
            ) + labs(title = paste(diff_cancers[j]))
          
          nes_plot_data <-
            nes_plot_data$data # extract the data from the plotEnrichment results
          
          nes_plot_data$cancer <-
            diff_cancers[j] # annotate the gene of interest
          
          # now, add the NES score calculated earlier
          nes <-
            cancer_enrichment$fgseaResult_cancers_depmap_final |>
            subset(pathway == diff_cancers[j]) |>
            select(NES) |>
            unique() |>
            as.numeric() |>
            round(digits = 2)
          
          nes_plot_data$NES <- nes
          
          # now, add the padj score calculated earlier
          padj <-
            cancer_enrichment$fgseaResult_cancers_depmap_final |>
            subset(pathway == diff_cancers[j]) |>
            select(padj) |>
            unique() |>
            as.numeric() |>
            round(digits = 2)
          
          nes_plot_data$padj <- padj
          
          # Bind the extracted data to the combined data frame
          combined_data2 <- rbind(combined_data2, nes_plot_data)
        }
        
        combined_data2$x <- as.numeric(combined_data2$x)
        combined_data2$y <- as.numeric(combined_data2$y)
        
        combined_data2 <- combined_data2 |>
          arrange(-NES) |>
          group_by(cancer) |>
          mutate(NES = replace(NES, row_number() > 1, NA)) |>
          ungroup()
        
        # Convert facet variable to factor with desired order
        combined_data2$cancer <-
          factor(combined_data2$cancer, levels = unique(combined_data2$cancer))
        
        # NES curves
        p10 <- ggplot(combined_data2,
                      aes(color = x,
                          x = x,
                          y = y)) +
          geom_hline(yintercept = 0) +
          geom_point(size = 2) +
          geom_line() +
          # geom_bar(aes(fill = y), stat = "identity", width = 0.01) +
          scale_color_viridis(name = "Cancer cell line\n gene effect rank", direction = 1) +
          scale_fill_viridis(direction = -1) +
          theme(
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = 'white'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title.align = 0.5,
            strip.text = element_text(size = 15)
          ) +
          xlab("Cell line rank order\n<- more dependent lines | less dependent lines ->") +
          ylab("Enrichment Score") +
          ggtitle(
            paste(
              "Most differentially enriched cancers by",
              genes_of_interest[i],
              "effect score"
            )
          ) +
          theme(plot.title = element_text(hjust = 0.5)) +
          guides(fill = FALSE) +
          facet_wrap( ~ cancer)
        
        # add NES scores to each plot
        nes_data <- combined_data2 |>
          subset(!is.na(NES))
        
          x_pos2 <- max(combined_data2$x) * 0.8
          y_pos2 <- max(combined_data2$y) * 0.8
        
        p10 <- p10 +
          geom_text(
            color = "black",
            size = 3.5,
            data = nes_data,
            mapping = aes(
              x = x_pos2,
              y = y_pos2,
              label = paste("NES = ", NES, "\n", "padj = ", padj, sep = "")
            )
          )
        
        # save data
        results_list[[9]] <- combined_data2
        
        # add plot to plot list
        plot_list[[10]] <- p10
      }
    }
    return(list(plots = plot_list, data = results_list))
  }