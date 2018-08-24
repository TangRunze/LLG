source("main_functions.R")
source("parse_neurodata.R")

# This function goes through a  set of parameters and computes
# the estimator and error
err_all_params <- function(alist_df, all_param, P){
    n <- nrow(alist_df$adj[[1]])
    abar <- Reduce("+", alist_df$adj) / nrow(alist_df)


    est_df <- purrrlyr::invoke_rows(compute_phat,
            all_param,
            alist = alist_df$adj,
            abar = abar,
            .to = "phat") %>%
        rowwise() %>%
        mutate(d = phat$dim, err = norm(phat$phat - P)^2 / n^2) %>%
        select(-phat) %>%
        bind_rows(tibble(d = n, dim = "abar", err = norm(abar - P)^2 / n^2))
    # group_by_at(vars(-phat)) %>%
    #   mutate(d = phat[[1]]$dim, phat = phat[[1]][1])
    #%>%
    #   mutate(err = phat %>% map_dbl(function(est) norm(est - P) ^ 2)) %>%
    #   select(-phat)
}

all_param <-  list(dim = c("ZG", "USVT"),
    diag_aug = c(TRUE,FALSE), diag_aug_sec = c(TRUE,FALSE),
    threshold = c(TRUE,FALSE), is_svd = c(TRUE,FALSE)) %>%
  cross_df

data_dir <- "~/Dropbox/Data/neurodata_dtmri/graphml/"

# Load data
data_df <- neurodata_df %>%
    filter(atlas == "CPAC200", dataset == "SWU4") %>%
    mutate(atlas = as.character(atlas)) %>%
    bind_cols(pmap_df(., load_graphml, data_dir = data_dir)) %>%
    group_by(dataset, atlas, subject, scan) %>%
    filter(is(adj[[1]], "Matrix")) %>%
    rowwise() %>% mutate(adj = list(adj > 0))

# compute P = average of all
P <- Reduce("+", data_df$adj) / nrow(data_df)

# 100 replicates
nmc <- 10
# for each m in 1,5,10, 20, 50
mrange <- cross_df(list(mc = 1:nmc, m = c(1, 2, 5, 10, 20)))

# all_param is in main_functions.R
# we get rid unthreshold because not thresholding is dumb
all_param <- all_param %>% filter(threshold)

mc <- 0
err_sample_df <-  mrange %>% group_by(mc) %>%
    mutate(res = map(m, function(m){
        print(mc[1])
        # sample m graphs
        sample_df <- data_df %>% ungroup() %>% sample_n(m)
        # compute all the errors
        err_df <- try(sample_df %>% err_all_params(all_param, P))
        # get rid of the adjacencies and store everything together
        sample_df %>% select(-adj) %>%
            nest(-atlas, -size, .key = "sample") %>%
            mutate(err = list(err_df))
    }))
 save(err_sample_df, file = "../../Data/test_all_param_desikan.RData")

load("../../Data/test_all_param_desikan.RData")
# Get the stuff we want
err_sample_df <- err_sample_df %>% unnest(res) %>% unnest(err)
err_sample_df %>% group_by(mc, m) %>% mutate(re = err / err[dim == "abar"])

# Make a few plots
abar_err_df <- err_sample_df %>% filter(dim == "abar") %>%
    ungroup %>% select(m, err)
phat_err_df <- err_sample_df %>% group_by(mc, m) %>%
    mutate(re = err / err[dim == "abar"]) %>% filter(dim != "abar")
phat_err_df %>% filter(threshold) %>%
    mutate(order = ifelse(is_svd, "LM", "LA")) %>%
    ggplot(aes(x = m, y = err, color = diag_aug, linetype = diag_aug_sec)) +
    stat_summary(fun.y = "mean", geom = "line") +
    stat_summary(data = abar_err_df,
        fun.y = "mean", geom = "line", color = "black", linetype = 1) +
    facet_grid(dim~order)

phat_err_df %>% filter(threshold) %>%
    mutate(order = ifelse(is_svd, "LM", "LA")) %>%
    ggplot(aes(x = m, y = re,
        color = diag_aug, linetype = diag_aug_sec,
        shape = diag_aug_sec)) +
    stat_summary(fun.y = "mean", geom="line") +
    stat_summary(fun.data = "mean_cl_boot") +
    geom_hline(yintercept = 1, color = "black", alpha = .3) +
    facet_grid(dim~order, scales = "free") +
    coord_cartesian(ylim = c(0.5, 2.5))

ggsave("../../ieee_tmi_response/compare_param_Desikan_SWU4.pdf",
    width = 9, height = 9)

phat_err_df %>% filter(diag_aug, diag_aug_sec) %>%
    ggplot(aes(x = d)) + geom_histogram() +
    facet_grid(m~is_svd + dim)

phat_err_df %>% filter(threshold, is_svd) %>%
    mutate(order = ifelse(is_svd, "LM", "LA")) %>%
    mutate(diag = diag_aug + 2 * diag_aug_sec,
        `Diag. Aug.` = recode(diag,
            `0` = "__", `1` = "_S", `2` =  "M_", `3` = "MS")) %>%
    ggplot(aes(x = d, y = re,
        linetype = `Diag. Aug.`, group = m + 400 * diag,
        color = m)) + #, linetype = diag_aug_sec)) +
    # stat_summary(fun.y = "mean", geom="line") +
    stat_smooth()  +
    geom_hline(yintercept = 1, color = "black", alpha = .3) +
    facet_grid(order ~ dim, scales = "free") +
    theme(legend.position = "bottom")

ggsave("../../ieee_tmi_response/compare_param_multi.pdf")

phat_err_df %>% filter(is_svd) %>%
    ggplot(aes(x = d, y = re, group = m,
        color = m)) +
    stat_smooth()  +
    geom_hline(yintercept = 1, color = "black", alpha = .3) +
    facet_grid(is_svd ~ dim, scales = "free") +
    theme(legend.position = "bottom")

phat_err_df %>% filter(threshold) %>%
    mutate(order = ifelse(is_svd, "LM", "LA")) %>%
    ggplot(aes(x = paste(dim, order, diag_aug, diag_aug_sec), y = re)) +
    stat_summary(fun.y = "mean", geom="point") +
    facet_grid(m~dim, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


# compute P
# compute all 32 possible
# also combute barA
# compute F-norm error to P