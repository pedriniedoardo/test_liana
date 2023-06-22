# Context Factorisation with tensor-cell2cell -----------------------------
# Daniel Dimitrov
# Saezlab, Heidelberg University
# daniel.dimitrov@uni-heidelberg.de
# 2023-02-24
# Source: vignettes/liana_cc2tensor.Rmd

# Introduction ------------------------------------------------------------
# Tensor decomposition as proposed in the tensor_cell2cell paper, enables us to decipher context-driven intercellular communication by simultaneously accounting for an unlimited number of “contexts”. These contexts could represent samples coming from longtidinal sampling points, multiple conditions, or cellular niches.

# The power of tensor-cell2cell is in its ability to decompose latent patterns of intercellular communication in an untargeted manner, in theory being able to handle cell-cell communication results coming from any experimental design, regardless of its complexity.

# Simply put, tensor_cell2cell uses LIANA’s output by sample to build a 4D tensor, represented by 1) contexts, 2) interactions, 3) sender, and 4) receiver cell types. This tensor is then decomposed into a set of factors, which can be interpreted as low dimensionality latent variables (vectors) that capture the CCC patterns across contexts.

# Here, we will combine LIANA with tensor_cell2cell to decipher potential ligand-receptor interaction changes. As a simple example, we will look at ~14000 PBMCs from 8 donors, each before and after IFN-beta stimulation (GSE96583; obtained via ExperimentHub & muscData). Note that by focusing on PBMCs, for the purpose of this tutorial, we assume that coordinated events occur among them.

# This tutorial was heavily influenced by the tutorials of tensor_cell2cell.

# Any usage of liana x tensor_cell2cell should logically cite both articles, and in particular tensor_cell2cell (see reference at the bottom).


# Load required libraries -------------------------------------------------
library(tidyverse)
library(SingleCellExperiment)
library(reticulate)
library(magrittr)
library(liana)
library(ExperimentHub)
library(patchwork)

# Request Data and Preprocess ---------------------------------------------

# Request Data ------------------------------------------------------------
# eh <- ExperimentHub()
# # Get Data
# sce <- eh[["EH2259"]]
# # save the object for future needs
# saveRDS(sce,"../../data/EH2259.rds")
sce <- readRDS("../../data/EH2259.rds")

# Preprocess --------------------------------------------------------------
# basic feature filtering
sce <- sce[rowSums(counts(sce) >= 1) >= 5, ]

# basic outlier filtering
qc <- scater::perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- scater::isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# Remove doublets
sce <- sce[, sce$multiplets=="singlet"]

# Set rownames to symbols
rownames(sce) <- rowData(sce)$SYMBOL

# log-transform
sce <- scuttle::logNormCounts(sce)

# Create a label unique for every sample
sce$context <- paste(sce$stim, sce$ind, sep="|")


# Ensure Consitency across Cell identities --------------------------------
# To obtain consistent CCC patterns across samples, we need to make sure that the cell identities are stable. Meaning there are enouhg cells per sample

# Plot
sce %>%
  get_abundance_summary(sample_col = "context",
                        idents_col = "cell", 
                        min_cells = 10, # min cells per sample
                        min_samples = 3, # min samples
                        min_prop = 0.2 # min prop of samples
  ) %>%
  plot_abundance_summary()

# filter non abundant celltypes
sce <- liana::filter_nonabundant_celltypes(sce,
                                           sample_col = "context",
                                           idents_col = "cell")

# replot the data
sce %>%
  get_abundance_summary(sample_col = "context",
                        idents_col = "cell", 
                        min_cells = 10, # min cells per sample
                        min_samples = 3, # min samples
                        min_prop = 0.2 # min prop of samples
  ) %>%
  plot_abundance_summary()

# Run liana for on each individual sample. --------------------------------
# In order to construct the tensor, we first need to obtain CCC predictions for each sample. In this case, we use SingleCellSignalR scores, as they are regularized, and theory directly comparable between any dataset. One can use any method with non-negative scores from LIANA /w cell2cell_tensor as they were previously shown to yield consistent results (Armingol & Baghdassarian, 2022).

# Note that liana_bysample works with SingleCellExperiment alone, if you wish to use Seurat, please use the as.SingleCellExperiment function.

# Run LIANA by sample
sce <- liana_bysample(sce = sce,
                      sample_col = "context",
                      idents_col = "cell",
                      method = "sca", # we use SingleCellSignalR's score alone
                      expr_prop = 0, # expression proportion threshold
                      inplace=TRUE, # saves inplace to sce
                      return_all = FALSE # whether to return non-expressed interactions 
)

# Note that in this case, we don’t apply the expr_prop threshold and we keep the values of all interactions as they are. Alternatively, one could be more restrictive with the communication events thought to be occurring, and increase expr_prop to a higher value, and return_all interactions with those below expr_prop being assigned to the worst inferred value. That value can also be manually replaced with e.g 0.

# We expect to see a successful LIANA run for each sample/context.
summary(sce@metadata$liana_res)

# Cell-cell Communication Tensor Decomposition ----------------------------
# Here, we call tensor_cell2cell.
# This function will first format the ligand-receptor scores per sample into a 4 Dimensional tensor.
# It will then estimate the number of factors to which the tensor will be decomposed (set rank to NULL, for the sake of computational speed, I have pre-calculated the rank and I explicitly set it to 7 here). Optimal rank estimation can be computationally demanding, but required to ensure robust results, one could consider increasing the runs parameter. See below an elbow plot example, calculated if rank is set to NULL.
# One can think of this as a higher-order non-negative matrix factorization where the factors can be used to reconstruct the tensor. We refer the user to the publication of tensor-cell2cell for further information.

# Note that by default LIANA will set up a conda environment with basilisk (if conda_env is NULL), the user can alternatively specify the name of a conda_env with cell2cell to be called via reticulate.

sce <- liana_tensor_c2c(sce = sce,
                        score_col = 'LRscore',
                        rank = 7,  # set to None to estimate for you data!
                        how='outer',  #  defines how the tensor is built
                        conda_env = "env_cell2cell", # used to pass an existing conda env with cell2cell
                        use_available = T # detect & load cell2cell if available
)

# The how parameter plays is essential in how we treat missing cell types and interactions. In this scenario, we use outer, which will only decompose CCC across any cell identities and interactions in the dataset. Alternatively, one could change it to ‘inner’ then CCC paterns will be decomposed patterns across cell identities & interactions shared between all samples/contexts. Similarly, one could consider how to fill missing values across the samples, namely missing interactions cell types via the lr_fill and cell_fill parameters, respecitvely. These are set to NaN by default and are thus imputed. For example, setting them to 0 would suggest that missing values are biologically-relevant.

# cell2cell_tensor also accepts additional parameters to further tune tensor decomposition, but these are out of scope for this tutorial. Stay tuned for a set of tutorials with comprehensive application of liana with tensor!

# Rank Estimation Example -------------------------------------------------
# Upon Rank estimation, an optimal rank will be estimated and returned according to the elbow of the normalized reconstruction error.

# This error measures how different is the original tensor from the sum of R rank-1 tensors used to approximate it. We want to pick a rank that provides a good trade-off that minimizes the number of factors (ranks) and the error.

# load pre-computed results 
# these would be assigned to this element, if computed
sce@metadata$tensor_res$elbow_metric_raw <- 
  read.csv(file.path(system.file(package = "liana"), "example_elbow.csv"))

# Estimate standard error
error_average <- sce@metadata$tensor_res$elbow_metric_raw %>%
  t() %>%
  as.data.frame() %>%
  mutate(rank=row_number()) %>% 
  pivot_longer(-rank, names_to = "run_no", values_to = "error") %>%
  group_by(rank) %>%
  summarize(average = mean(error),
            N = n(),
            SE.low = average - (sd(error)/sqrt(N)),
            SE.high = average + (sd(error)/sqrt(N))
  )

# plot
error_average %>%
  ggplot(aes(x=rank, y=average), group=1) +
  geom_line(col='red') + 
  geom_ribbon(aes(ymin = SE.low, ymax = SE.high), alpha = 0.1) +
  geom_vline(xintercept = 7, colour='darkblue') + # rank of interest
  theme_bw() +
  labs(y="Error", x="Rank")

# Results Overview --------------------------------------------------------
# Upon a seccusseful run cell2cell_tensor will decompose the tensor into a set of factors, or four vectors corresponding to the initial dimensions of the tensor.
# contexts - the factor scores assigned to each sample/context
# interactions - the interaction loadings per factor
# senders - loadings for sender/source cell identities
# receivers - loadings for receivers/target cell identities

# get the factors
factors <- get_c2c_factors(sce, group_col = "stim", sample_col = "context")

# show them
glimpse(factors)

# Here, we examine the behavior of the different dimensions across the factors. When we look at the contexts (samples) loadings Factor.4 seems to be notably different between the groups. We also see that the loadings for the “Sender” and “Receiver” cells have a relatively uniform distribution, suggesting that most cell types are involved in the CCC events that distinguish the conditions.

# Plot overview
plot_c2c_overview(sce, group_col="stim", sample_col="context")

# Statistical comparison of Communication Patterns ------------------------
# Here, we see that the communication patterns (or context loadings) identified statistically significant patterns before and after stimulation.

# These factors thus represent differences in the ligand-receptor interactions as well as cell types participating in cell-cell communication before and after IFN-beta stimulation.

# Get all boxplots
all_facts_boxes <- plot_context_boxplot(sce,
                                        sample_col = "context",
                                        group_col = "stim",
                                        test="t.test", # applicable only to two groups
                                        paired=TRUE #! Is this the case for your data?
)

# Combine all boxplots
# require(patchwork)
wrap_plots(
  all_facts_boxes,
  ncol=4) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom') +
  theme(plot.tag = element_text(face = 'bold', size = 16)
  )

# Contexts Heatmap --------------------------------------------------------
# Again, we see a clear seperation between STIM and CTRL, further suggesting the relevance of the changes in the inferred interactions.
plot_context_heat(sce, group_col = "stim", sample_col="context")

# Cell-cell pairs with high potential of interaction ----------------------
# Here we use the product of the source and target loadings to generate the heatmap of potential celltype-celltype pair relationships, which contribute to Factor.4.
plot_c2c_cells(sce,
               factor_of_int = "Factor.4",
               name = "Loadings \nProduct")

# Gini Coefficintes -------------------------------------------------------
# Gini Coefficients of Factor-specific Communicating Sender and Receiver cell types
# Gini coefficients range from 0 to 1, and are measures of dispersion, typically used to measure inequality.

# Here, the Gini coefficient is used to measure the imbalance of communication in terms of the cell types.

# Gini coefficient of 1 would suggests that there is a single cell type contributing to the communication patterns within a factor, while a value of 0 suggest perfect equality between the cell type.

# When we focus on Factor.4, we see that the both source/sender and target/receiver cell types have nearly uniform contributions.

# Get loadings for source/sender Cell types
calculate_gini(factors$senders)

# Get loadings for target/receiver Cell types
calculate_gini(factors$receivers)

# LR loadings Heatmap -----------------------------------------------------
# We can see the LRs involved across contexts.

# Though, perhaps since in this case Factors 4 and 5 are those associated with the stimulation, perhaps it’s best to focus on those.
plot_lr_heatmap(sce,  n = 5, cluster_columns=FALSE)

# LR loadings Pathway Enrichment ------------------------------------------
# Footprint Enrichment Analysis

# Here, we will use pathways from PROGENy, together with decoupleR, to do an enrichment analysis on the LR loadings.

# PROGENy represent a data-driven signatures for 14 pathways, while decoupleR is an framework with multiple enrichment methods.

# Since PROGENy was originally intended to be used with genes, we need to reformat it to match ligand-receptor predictions. Namely, generate_lr_geneset will assign a LR to a specific pathway, only if both the ligand and receptor, as well as all their subunits, are associated with the same pathway. And also in the case of weighted gene sets (like progeny), all entities have the same weight sign.

# obtain progeny gene sets
progeny <- decoupleR::get_progeny(organism = 'human', top=5000) %>%
  select(-p_value)

# convert to LR sets
progeny_lr <- generate_lr_geneset(sce,
                                  resource = progeny)

progeny_lr

# Enrichment Dotplot ------------------------------------------------------
# Here, we see that LRs associated with the JAK-STAT gene in progeny are enriched in Factor.4.

# interaction loadings to matrix
mat <- factors$interactions %>%
  column_to_rownames("lr") %>%
  as.matrix()

# run enrichment analysis with decoupler
# (we fit a univariate linear model for each gene set)
# We don't consider genesets with minsize < 10
res <- decoupleR::run_ulm(mat = mat,
                          network = progeny_lr,
                          .source = "set",
                          .target = "lr",
                          minsize=10) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

res %>% # sig/isnig flag
  mutate(significant = if_else(p_adj <= 0.05, "signif.", "not")) %>%
  ggplot(aes(x=source, y=condition, shape=significant,
             colour=score, size=-log10(p_value+1e-36))) +
  geom_point() +
  scale_colour_gradient2(high = "red", low="blue") +
  scale_size_continuous(range = c(3, 12)) +
  scale_shape_manual(values=c(21, 16)) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Pathway",
       y="Factor",
       colour="Activity"
  )


# LRs driving enrichment of JAK-STAT in Factor.4 --------------------------
# Plot LRs associated with Estrogen
lrs <-  factors$interactions %>%
  left_join(progeny_lr, by="lr") %>%
  filter(set=="JAK-STAT") %>%
  select(lr, set, mor, loading = Factor.4) %>%
  mutate(lr = gsub(as.character(str_glue("\\^")), " -> ", lr)) %>%
  mutate(weight = if_else(mor >= 0, "positive", "negative"))
lrs %>%
  # only label those that are > x
  mutate(lr = if_else(loading>=0.001 & abs(mor) > 2, lr, "")) %>%
  ggplot(aes(x=mor, y=loading, colour=weight)) +
  # label only top 20
  stat_smooth(method = "lm", col = "red") +
  geom_point(alpha = 0.5) + 
  ggrepel::geom_label_repel(aes(label = lr)) +
  theme_bw(base_size = 15) +
  scale_colour_manual(values = c("royalblue3", "red")) +
  labs(x="Pathway Weight", y="LR Loading")
