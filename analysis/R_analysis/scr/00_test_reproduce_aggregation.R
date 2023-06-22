# basic functions ---------------------------------------------------------
# source 
# https://github.com/saezlab/liana/blob/86adaca8832bc2beb572d4d83d837e480d967865/R/liana_aggregate.R#L196
.rank_matrix <- function(glist, ...){
  
  # Get unique entities
  u.ents <- unique(c(glist, recursive = TRUE))
  
  # num of entities per col
  num.ents <- length(u.ents)
  
  # get position for each vect
  pos.mat <- sapply(FUN = match,
                    X = glist,
                    x = u.ents,
                    nomatch = num.ents) %>%
    matrix(nrow = num.ents,
           ncol = length(glist),
           dimnames = list(u.ents, names(glist)))
  
  # Fill mat /w total ents
  rank.mat <- matrix(num.ents,
                     nrow = num.ents,
                     ncol = length(glist),
                     dimnames = list(u.ents, names(glist)))
  
  # return rank/total_rank by pos
  return(pos.mat / rank.mat)
  
}

.robust_rank_agg <- function(rmat){
  tibble(interaction = rownames(rmat), # Name
         # calc aggr ranks
         aggregate_rank = unname(apply(rmat, 1, .rho_scores))) %>% # Score
    arrange(aggregate_rank)
}

.rho_scores <- function(r){
  r <- sort(r)
  n <- length(r) #length is sometimes larger than max -> over-inflates FPs
  
  # Calc beta p-vals
  p <- pbeta(q=r,
             shape1 = 1:n,
             shape2 = n - (1:n) + 1)
  
  # correct beta pvals
  .corr_beta_pvals(p = min(p), k=n)
}

.corr_beta_pvals <- function(p, k){
  min(p * k, 1)
}


# -------------------------------------------------------------------------
# main function
.aggregate_rank <- function(liana_mlist, join_cols, verbose, ...){
  # liana_message("Aggregating Ranks", output = "message", verbose = verbose)
  
  liana_mlist %>%
    map(function(res){
      # bad practice, but almost unavoidable here...
      res %>%
        unite(c("source", "target",
                "ligand.complex", "receptor.complex"),
              col = "interaction",
              sep = "⊎") %>%
        pull("interaction")
    }) %>%
    .rank_matrix %>%
    .robust_rank_agg(.,
                     ...) %>%
    separate(col = "interaction", sep = "⊎",
             into = c("source", "target",
                      "ligand.complex", "receptor.complex"))
}

# -------------------------------------------------------------------------
# build the test set
liana_path <- system.file(package = "liana")
testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
# run liana
liana_res <- liana_wrap(testdata, method=c("natmi","sca","logfc"))

# build the final result to confirm the approach
test_rank <- liana_res %>% liana_aggregate()

# extract only the ranks from the final result
list_complex_test <- lapply(names(liana_res),function(x){
  test_rank %>% 
    select(source,target,ligand.complex,receptor.complex,contains(x)) %>% 
    # arrange(pick(ends_with("rank")))
    arrange(pick(paste0(x,".rank")))
}) %>% set_names(names(liana_res))

# try to replicate the average rank
test_average <- list_complex_test %>% 
  purrr::reduce(left_join,by=c("source","target","ligand.complex","receptor.complex")) %>% 
  rowwise() %>% 
  mutate(avg = mean(c(natmi.rank, sca.rank, logfc.rank))) %>% 
  arrange(avg) %>% 
  select(source,target,ligand.complex,receptor.complex,avg)

# compare with the original result
left_join(test_rank,test_average,by = c("source","target","ligand.complex","receptor.complex")) %>% 
  mutate(test = mean_rank == avg) %>% 
  summarise(test = sum(!test))
# confirmed. the mean_rank is the average of the individula ranks
  
# test replicate the aggregation
test_aggregate <- .aggregate_rank(list_complex_test)

# compare it with the original reuslt
left_join(test_rank,test_aggregate,by = c("source","target","ligand.complex","receptor.complex")) %>% 
  mutate(test = aggregate_rank.x == aggregate_rank.y) %>% 
  summarise(test = sum(!test))
# sum(!test_rank$aggregate_rank == test_aggregate$aggregate_rank)
# confirmed. the aggregate_rank can be reproduced running the .aggregate_rank on the list object

# -------------------------------------------------------------------------
# try to understand the individual steps
test_01 <- list_complex_test %>%
  map(function(res){
    # bad practice, but almost unavoidable here...
    res %>%
      unite(c("source", "target",
              "ligand.complex", "receptor.complex"),
            col = "interaction",
            sep = "⊎") %>%
      pull("interaction")
  })

# the output of this step is the relative rank devided by the total number of rows of the dataaset (total number of pairs)
test_02 <- test_01 %>% 
  .rank_matrix

# glist <- test_01
# 
# # Get unique entities
# u.ents <- unique(c(glist, recursive = TRUE))
# 
# # num of entities per col
# num.ents <- length(u.ents)
# 
# # get position for each vect
# pos.mat <- sapply(FUN = match,
#                   X = glist,
#                   x = u.ents,
#                   nomatch = num.ents) %>%
#   matrix(nrow = num.ents,
#          ncol = length(glist),
#          dimnames = list(u.ents, names(glist)))
# 
# # Fill mat /w total ents
# rank.mat <- matrix(num.ents,
#                    nrow = num.ents,
#                    ncol = length(glist),
#                    dimnames = list(u.ents, names(glist)))
# 
# # return rank/total_rank by pos
# pos.mat / rank.mat

# final ranking
test_03 <- test_02 %>%
  .robust_rank_agg

data.frame(res = apply(test_02, 1, .rho_scores)) %>% 
  arrange(res)

# .rho_scores <- function(r){
#   r <- sort(r)
#   n <- length(r) #length is sometimes larger than max -> over-inflates FPs
#   
#   # Calc beta p-vals
#   p <- pbeta(q=r,
#              shape1 = 1:n,
#              shape2 = n - (1:n) + 1)
#   
#   # correct beta pvals
#   .corr_beta_pvals(p = min(p), k=n)
# }
# 
# .corr_beta_pvals <- function(p, k){
#   min(p * k, 1)
# }
