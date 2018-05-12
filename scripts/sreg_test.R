library(pbr)
library(ggplot2)

soy <- gauch.soy

sreg_rss <- combn(x = levels(soy$env), m = 2, function(env_pair) {
  fit <- sreg(yield ~ gen + env + gen:env, data = subset(soy, env %in% env_pair))
  rss <- sum(fit$effects$residual^2)
  data.frame(env1 = env_pair[1], env2 = env_pair[2], rss = rss, row.names = NULL, stringsAsFactors = FALSE)
}, simplify = FALSE)

# Fit using all data
sreg_all <- sreg(yield ~ gen + env + gen:env, data = soy)

# Combine the data.frame
sreg_rss_df <- sreg_rss %>%
  bind_rows() %>%
  arrange(env1, env2)


# Convert to distance matrix
rss_mat <- xtabs(rss ~ env1 + env2, sreg_rss_df)
rss_dist <- as.dist(t(rss_mat))


## Plot all



## Clusters should group the environments according to the names in this vector:
env_group <- sreg_all$effects %>% distinct(env, gen, fitted) %>% group_by(env) %>% top_n(n = 1) %>% {setNames(nm = .$env, object = .$gen)}



rss_clust <- hclust(d = rss_dist, method = "complete")

sreg_cut <- cutree(tree = rss_clust, k = 2)


# Combine
env_group_df <- merge(x = as.data.frame(sreg_cut), y = as.data.frame(env_group), by = "row.names")


