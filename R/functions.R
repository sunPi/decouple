#---- Functions ----
plot_pathways       <- function(df, test){
  colors <- rev(RColorBrewer::brewer.pal(n = nrow(df), name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = df, 
                       mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                              y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 10, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    ggplot2::labs(x = "Pathways", y = "Enrichment Score", title = test) +
    coord_flip()+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}
plot_single_pathway <- function(net, deg, pathway){
  
  df <- net %>%
    dplyr::filter(source == pathway) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(ID = target, 
                  color = "3") %>%
    tibble::column_to_rownames('target')
  
  inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
  
  df <- df[inter, ]
  
  df['t_value'] <- deg[inter, ]
  
  df <- df %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = df, 
                       mapping = ggplot2::aes(x = weight, 
                                              y = t_value, 
                                              color = color)) + 
    ggplot2::geom_point(size = 2.5, 
                        color = "black") + 
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
    ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
    ggplot2::ggtitle(pathway)
  
  return(p)
}
getColData          <- function(exp.data, cts, reference){
  coldata <- exp.data[order(exp.data[[1]]), ]
  names(coldata)[2] <- "condition"
  coldata$condition <- relevel(factor(coldata$condition), ref = reference)
  
  return(coldata)
}
preprocess <- function(cts, exp, reference){
  # Check if in ENSEMBL format and convert into HUGO
  if(all(grepl(pattern = "^ENSG*", cts$gene_id))){
    # gene.list$gene_id <- mapIds(org.Hs.eg.db, keys = row.names(gene.list), keytype = "ENSEMBL", column = "ENTREZID")
    cts$gene_id <- mapIds(org.Hs.eg.db, keys = cts$gene_id, keytype = "ENSEMBL", column = "SYMBOL")
  }
  
  # Remove duplicates and add row names if you want to normalize counts without deseq use "cts <- log2(cts +1)"
  if(any(duplicated(cts$gene_id))){
    print("Found duplicated genes, removing them...")
    cts <- cts[-which(duplicated(cts$gene_id)),]
  }
  
  row.names(cts) <- NULL
  gene.id <- cts[,1]
  cts     <- as.matrix(cts[,-1])
  row.names(cts) <- gene.id
  
  col <- getColData(exp, cts, reference = reference)
  
  # setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
  setClass("ExpData", contains = "SummarizedExperiment")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                        colData = col,
                                        design = ~ condition)
  
  dds <- DESeq2::estimateSizeFactors(dds)
  norm.counts <- DESeq2::counts(dds, normalized = TRUE)
  cts.log2.norm <- log2(norm.counts + 1)
  
  if(any(duplicated(deseq.res$deg$gene_id))){
    print("Found duplicated genes, removing them...")
    deseq.res$deg <- deseq.res$deg[-which(duplicated(deseq.res$deg$gene_id)),]
  }
  
  # MesoLab colour Scheme
  # #5f0053 - Downregulated
  # #087c95 - Upregulated
  
  # Create a decoupleR deg dataframe
  deg <- data.frame(row.names = deseq.res$deg$gene_id, "stat"= deseq.res$deg$stat)
  
  return(list(normalized.counts = cts.log2.norm,
              deg               = deg))
}
split_decouple <- function(results, outfolder, do.split){
  if(do.split){
    # Split the dataframe based on the 'statistic' column
    split.res <- split(results, results$statistic)
    
    # Remove repeated values in each data frame of a list
    split.res <- lapply(split.res, function(df) df[, -c(1:3)])
    
    for(n in names(split.res)){
      df <- split.res[[n]]
      out <- here(outfolder, paste0(n,".csv"))
      
      # out <- here(outfolder)
      write.csv(df, out, row.names = F)
    }
    
  } else{
    
    results <- results[, -c(1:2)]
    results <- results[results$p_value < 0.05,]
    
    out <- here(outfolder, paste0(methods,".csv"))
    # out <- here(outfolder)
    write.csv(results, out, row.names = F)
    
  }
}
get_consensus <- function(results, outfolder){
  
  res.consensus <- run_consensus(results)
  
  if(any(res.consensus$p_value < 0.05)){
    res.sig <- res.consensus[res.consensus$p_value < 0.05,]
    res.sig <- res.sig[order(res.sig$source),]
    
    write.csv(res.sig, here(outfolder, "consensus.csv"))
    return(res.sig)
    
  } else{
    message("There are no significant consensus scores!")
    
  }
}
decouple_pathways <- function(cts, deg, top.n, methods, do.split, outfolder){
  message("Running activity inferrence using the PROGENy with ", top.n, " significant associations...")
  net <- decoupleR::get_progeny(organism = 'human',
                                top = top.n)
  
  # Decouple Pathways
  results <- decouple(
    mat     = cts,
    network = net,
    .source = "source",
    .target = "target",
    statistics = methods,
    consensus_score = F,
    minsize = 5
  )
  
  out.pth <- here(outfolder, "pathways")
  dir.create(out.pth, T, T)
  
  # Process and Save Decoupled Results
  split_decouple(results, out.pth, do.split)
  
  # Compute Consensus score
  message("Computing consensus of all methods...")
  consensus <- get_consensus(results, out.pth)
  
  message("Running contrast pathway analysis using MLM...")
  # Run mlm
  contrast_acts <- decoupleR::run_mlm(mat     = deg, 
                                      net     = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor    = 'weight', 
                                      minsize = 5)
  
  test <- unique(exp[!exp[,2] %in% reference,][,2])
  
  # if(methods == "fgsea"){
  #   results <- dplyr::filter(results, results$statistic %in% "norm_fgsea")
  #   contrast_acts <- dplyr::filter(contrast_acts, contrast_acts$statistic %in% "norm_fgsea")
  #   contrast.sig <- contrast_acts[contrast_acts$p_value < 0.05,]
  #   
  #   contrast.sig <- contrast.sig %>% 
  #     dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
  #   
  # } else{
  #   contrast.sig <- contrast_acts[contrast_acts$p_value < 0.05,]
  #   contrast.sig <- contrast.sig %>% 
  #     dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
  # }
  # 
  
  contrast.sig <- contrast_acts[contrast_acts$p_value < 0.05,]
  contrast.sig <- contrast.sig %>% 
    dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
  
  # Plot
  pathways <- plot_pathways(contrast.sig, test)
  
  message("Saving plots...")
  ggsave(here(out.pth, "pathways_waterfall.jpg"), plot = pathways,
         width = 8, height = 6, units = "in",
         dpi = 600, device = "jpeg")
  message("Finished!")
}
decouple_TF <- function(cts, deg, methods, do.split, outfolder){
  message("Running TF inferrence using the CollecTRI...")
  net <- decoupleR::get_collectri(organism = 'human', 
                                  split_complexes = FALSE)
  # Run decouple
  results <- decouple(
    mat = cts,
    network = net,
    .source = "source",
    .target = "target",
    statistics = methods,
    consensus_score = F,
    minsize = 5
  )
  
  out.tf  <- here(outfolder, "TF")
  dir.create(out.tf, T, T)
  
  # Process and Save Decoupled Results
  split_decouple(results, out.tf, do.split)
  
  # Compute Consensus score
  message("Computing consensus of all methods...")
  consensus <- get_consensus(results, out.tf)
  
  message("Running contrast TF analysis using ULM...")
  # Run ulm
  sample_acts <- decoupleR::run_ulm(mat = cts, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'mor', 
                                    minsize = 5)
  
  # TF Inferrence between contrasts
  n_tfs <- 25
  # Run ulm
  contrast_acts <- decoupleR::run_ulm(mat = deg[, 'stat', drop = FALSE], 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor='mor', 
                                      minsize = 5)
  
  # Filter top TFs in both signs
  f_contrast_acts <- contrast_acts %>%
    dplyr::mutate(rnk = NA)
  
  msk <- f_contrast_acts$score > 0
  
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  
  tfs <- f_contrast_acts %>%
    dplyr::arrange(rnk) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  f_contrast_acts <- f_contrast_acts %>%
    dplyr::filter(source %in% tfs)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  tf.p <- ggplot2::ggplot(data = f_contrast_acts, 
                          mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                                 y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 20),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 16, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 16, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    coord_flip() +
    ggplot2::labs(x = "Transcription", y = "Enrichment", title = test)
  
  message("Saving plots...")
  ggsave(here(out.tf, "tf_waterfall.jpg"), plot = tf.p,
         width = 8, height = 6, units = "in",
         dpi = 300, device = "jpeg")
  message("Finished!")
}
# Run decouple
# RunDC <- function(cts, top.n, methods, do.split, outfolder){
#   
#   
#   
#   
#   ### Run Decouple on PROGENy (Pathway RespOnsive GENes for Activity Inference)
#   decouple_pathways(cts, top.n, methods, do.split, outfolder)
#   ### Run Decouple for TF Inference
#   decouple_TF(cts, top.n, methods, do.split, outfolder)
# }

# 
# 
# # results <- decouple(
# #   mat = cts.log2.norm,
# #   network = net,
# #   .source = "source",
# #   .target = "target",
# #   statistics = methods,
# #   consensus_score = F,
# #   minsize = 5
# # )
# # 
# # 
# # if(do.split){
# #   # Split the dataframe based on the 'statistic' column
# #   split.res <- split(results, results$statistic)
# #   
# #   # Remove repeated values in each data frame of a list
# #   split.res <- lapply(split.res, function(df) df[, -c(1:3)])
# #   
# #   for(n in names(split.res)){
# #     df <- split.res[[n]]
# #     out <- here(outfolder, paste0(n,".csv"))
# #     # out <- here(outfolder)
# #     write.csv(df, out, row.names = F)
# #   }
# #   
# # } else{
# #   results <- results[, -c(1:2)]
# #   results <- results[results$p_value < 0.05,]
# #   
# #   out <- here(outfolder, paste0(methods,".csv"))
# #   # out <- here(outfolder)
# #   write.csv(df, out, row.names = F)
# #   
# # }
# # 
# # # Compute Consensus score
# # message("Computing consensus of all methods...")
# # res.consensus <- run_consensus(results)
# # 
# # if(any(res.consensus$p_value < 0.05)){
# #   res.sig <- res.consensus[res.consensus$p_value < 0.05,]
# #   res.sig <- res.sig[order(res.sig$source),]
# #   
# #   write.csv(res.sig, here(outfolder, "consensus.csv"))
# #   
# # } else{
# #   message("There are no significant consensus scores!")
# # }
# 
# decoupled <- RunDC(cts.log2.norm, net, methods, do.split)
# 
# # Run mlm
# contrast_acts <- decoupleR::run_mlm(mat     = deg, 
#                                     net     = net, 
#                                     .source = 'source', 
#                                     .target = 'target',
#                                     .mor    = 'weight', 
#                                     minsize = 5)
# 
# test <- unique(exp[!exp[,2] %in% reference,][,2])
# if(methods == "fgsea"){
#   results <- dplyr::filter(results, results$statistic %in% "norm_fgsea")
#   contrast_acts <- dplyr::filter(contrast_acts, contrast_acts$statistic %in% "norm_fgsea")
#   contrast.sig <- contrast_acts[contrast_acts$p_value < 0.05,]
#   
#   contrast.sig <- contrast.sig %>% 
#     dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
#   
# } else{
#   contrast.sig <- contrast_acts[contrast_acts$p_value < 0.05,]
#   contrast.sig <- contrast.sig %>% 
#     dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
# }
# 
# # Plot
# pathways <- plot_pathways(contrast.sig, test)
# 
# ggsave(here(outfolder, "pathways_waterfall.jpg"), plot = pathways,
#        width = 8, height = 6, units = "in",
#        dpi = 600, device = "jpeg")
# 
# ### Run Decouple for TF Inference
# message("Running TF inferrence using the CollecTRI with...")
# net <- decoupleR::get_collectri(organism = 'human', 
#                                 split_complexes = FALSE)
# # Run decouple
# # results <- decouple(
# #   mat = cts.log2.norm,
# #   network = net,
# #   .source = "source",
# #   .target = "target",
# #   statistics = methods,
# #   consensus_score = F,
# #   minsize = 5
# # )
# 
# # Run ulm
# sample_acts <- decoupleR::run_ulm(mat = cts, 
#                                   net = net, 
#                                   .source = 'source', 
#                                   .target = 'target',
#                                   .mor = 'mor', 
#                                   minsize = 5)
# 
# 
# # TF Inferrence between contrasts
# n_tfs <- 25
# # Run ulm
# contrast_acts <- decoupleR::run_ulm(mat = deg[, 'stat', drop = FALSE], 
#                                     net = net, 
#                                     .source = 'source', 
#                                     .target = 'target',
#                                     .mor='mor', 
#                                     minsize = 5)
# 
# # Filter top TFs in both signs
# f_contrast_acts <- contrast_acts %>%
#   dplyr::mutate(rnk = NA)
# 
# msk <- f_contrast_acts$score > 0
# 
# f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
# f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
# 
# tfs <- f_contrast_acts %>%
#   dplyr::arrange(rnk) %>%
#   head(n_tfs) %>%
#   dplyr::pull(source)
# 
# f_contrast_acts <- f_contrast_acts %>%
#   dplyr::filter(source %in% tfs)
# 
# colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
# 
# tf.p <- ggplot2::ggplot(data = f_contrast_acts, 
#                      mapping = ggplot2::aes(x = stats::reorder(source, score), 
#                                             y = score)) + 
#   ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
#                     color = "black",
#                     stat = "identity") +
#   ggplot2::scale_fill_gradient2(low = colors[1], 
#                                 mid = "whitesmoke", 
#                                 high = colors[2], 
#                                 midpoint = 0) + 
#   ggplot2::theme_minimal() +
#   ggplot2::theme(axis.title = element_text(face = "bold", size = 20),
#                  axis.text.x = ggplot2::element_text(angle = 45, 
#                                                      hjust = 1, 
#                                                      size = 16, 
#                                                      face = "bold"),
#                  axis.text.y = ggplot2::element_text(size = 16, 
#                                                      face = "bold"),
#                  panel.grid.major = element_blank(), 
#                  panel.grid.minor = element_blank()) +
#   coord_flip() +
#   ggplot2::labs(x = "Transcription", y = "Enrichment", title = test)
# 
# ggsave(here(outfolder, "tf_waterfall.jpg"), plot = tf.p,
#        width = 8, height = 6, units = "in",
#        dpi = 300, device = "jpeg")
# 
# 
# 
# 
# 
# # Functions
# plot_wilcox <- function(df, x, y, fill, p.value, title, y.axis.label = NULL){
#   get_stars <- function(p) {
#     if (p < 0.001) return("***")
#     else if (p < 0.01) return("**")
#     else if (p < 0.05) return("*")
#     else return("ns")
#   }
#   p_star <- get_stars(p.value)
#   
#   # Plot boxplot
#   y_max <- max(y)
#   if(is.null(y.axis.label)){
#     y.name <- sub(".*\\$", "", deparse(substitute(y)))
#   } else{
#     y.name = y.axis.label
#   }
#   
#   
#   # groups[1]
#   # col.scale <- c("Q1" = "#077f97", "Q4" = "#400257")
#   # Plot with annotate
#   p <- ggplot(df, aes(x = x, y = y, fill = fill)) +
#     geom_boxplot(linewidth = 1.2) +
#     scale_fill_manual(values = c("#077f97", "#400257")) +
#     annotate("segment", x = 1, xend = 2, y = y_max + 0.2, yend = y_max + 0.2) +
#     annotate("text", x = 1.5, y = y_max + 0.4, label = p_star, size = 6) +
#     theme_prism(base_size = 20) +
#     labs(title=title,
#          x ="",
#          y = y.name)
#   return(p)
# }

# Code To Recycle
# Transform to wide matrix
# sample_acts_mat <- sample_acts %>%
#   tidyr::pivot_wider(id_cols = 'condition', 
#                      names_from = 'source',
#                      values_from = 'score') %>%
#   tibble::column_to_rownames('condition') %>%
#   as.matrix()
# 
# # Scale per feature
# sample_acts_mat <- scale(sample_acts_mat)
# sample_acts_mat
# # Color scale
# colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
# colors.use <- grDevices::colorRampPalette(colors = colors)(100)
# 
# my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
#                seq(0.05,2, length.out = floor(100 / 2)))
# 
# # Plot
# path.heat.p <- pheatmap::pheatmap(mat = sample_acts_mat,
#                    color = colors.use,
#                    border_color = "white",
#                    breaks = my_breaks,
#                    cellwidth = 20,
#                    cellheight = 20,
#                    treeheight_row = 20,
#                    treeheight_col = 20)
# 
# ggsave(here(outfolder, "pathways_heatmap.jpg"), plot = path.heat.p,
#        width = 8, height = 6, units = "in",
#        dpi = 600, device = "jpeg")
# plot_single_tf <- function(net, deg, tf){
#   df <- net %>%
#     dplyr::filter(source == tf) %>%
#     dplyr::arrange(target) %>%
#     dplyr::mutate(ID = target, color = "3") %>%
#     tibble::column_to_rownames('target')
#   
#   inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
#   
#   df <- df[inter, ]
#   
#   df[,c('logfc', 't_value', 'p_value')] <- deg[inter, ]
#   
#   df <- df %>%
#     dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
#     dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
#     dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
#     dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))
#   
#   colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
#   
#   p <- ggplot2::ggplot(data = df, 
#                        mapping = ggplot2::aes(x = logfc, 
#                                               y = -log10(p_value), 
#                                               color = color,
#                                               size = abs(mor))) + 
#     ggplot2::geom_point(size = 2.5, 
#                         color = "black") + 
#     ggplot2::geom_point(size = 1.5) +
#     ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
#     ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,
#                                                      size = 1)) + 
#     ggplot2::theme_minimal() +
#     ggplot2::theme(legend.position = "none") +
#     ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
#     ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
#     ggplot2::ggtitle(tf)
#   
#   return(p)
# }
# 
# foxi1.p <- plot_single_tf(net, deg, tf = "FOXI1")
# nfkb.p <- plot_single_tf(net, deg, tf = "NFKB")
# myc.p <- plot_single_tf(net, deg, tf = "MYC")
# # Transform to wide matrix
# sample_acts_mat <- sample_acts %>%
#   tidyr::pivot_wider(id_cols = 'condition', 
#                      names_from = 'source',
#                      values_from = 'score') %>%
#   tibble::column_to_rownames('condition') %>%
#   as.matrix()
# 
# # Get top tfs with more variable means across clusters
# tfs <- sample_acts %>%
#   dplyr::group_by(source) %>%
#   dplyr::summarise(std = stats::sd(score)) %>%
#   dplyr::arrange(-abs(std)) %>%
#   head(n_tfs) %>%
#   dplyr::pull(source)
# 
# sample_acts_mat <- sample_acts_mat[,tfs]
# 
# # Scale per sample
# sample_acts_mat <- scale(sample_acts_mat)
# 
# # Choose color palette
# colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
# colors.use <- grDevices::colorRampPalette(colors = colors)(100)
# 
# my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
#                seq(0.05, 2, length.out = floor(100 / 2)))
# 
# # Plot
# heat.tf.p <- pheatmap::pheatmap(mat = sample_acts_mat,
#                    color = colors.use,
#                    border_color = "white",
#                    breaks = my_breaks,
#                    cellwidth = 15,
#                    cellheight = 15,
#                    treeheight_row = 20,
#                    treeheight_col = 20)
# tgfb.p     <- plot_single_pathway(net, deg, pathway = "TGFb")
# estrogen.p <- plot_single_pathway(net, deg, pathway = "Estrogen")
# p53.p      <- plot_single_pathway(net, deg, pathway = "p53")


# ggsave(here(outfolder, "tgfb_coord.jpg"), plot = tgfb.p,
#        width = 8, height = 6, units = "in",
#        dpi = 600, device = "jpeg")
# ggsave(here(outfolder, "estrogen_coord"), plot = estrogen.p,
#        width = 8, height = 6, units = "in",
#        dpi = 600, device = "jpeg")
# ggsave(here(outfolder, "p53_coord"), plot = p53.p,
#        width = 8, height = 6, units = "in",
#        dpi = 600, device = "jpeg")