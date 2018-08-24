source("main_functions.R")
source("parse_neurodata.R")

data_dir <- "~/Dropbox/Data/neurodata_dtmri/graphml/"

aoi <- "CPAC200"
doi <- "SWU4"

# Load data
data_df <- neurodata_df %>%
    filter(atlas == aoi, dataset == doi) %>%
    mutate(atlas = as.character(atlas)) %>%
    bind_cols(pmap_df(., load_graphml, data_dir = data_dir)) %>%
    group_by(dataset, atlas, subject, scan) %>%
    filter(is(adj[[1]], "Matrix")) %>%
    rowwise() %>% mutate(adj = list(adj > 0))

# compute P = average of all
m_all <- nrow(data_df)
P <- Reduce("+", data_df$adj) / m_all

n <- nrow(P)
p_dhat_p <- compute_phat(data_df$adj)
d_p <- p_dhat_p$dim
p_dhat_p <- p_dhat_p$phat

p_dec <- decompose(P, n, is_svd = FALSE)

m <- 50
alist <- sample_n(data_df, m)$adj
abar <- Reduce("+", alist) / m



compare_to_latent <- function(alist, p, p_dec, p_dhat_p, d_p){
    m <- length(alist)
    abar <- Reduce("+", alist) / m

    # Compute estimate
    phat <- compute_phat(alist)
    dim_est <- phat$dim
    phat  <- phat$phat

    p_dhat <- compute_phat(alist, abar = p, dim = dim_est)$phat
    # low_rank_from_dec(p_dec, dim_est)

    # Compare Phat to P
    bind_rows(
        list(est = "abar", param = "p", de = n, dp = n,
            err = Matrix::norm(abar - p, "f")),
        list(est = "abar", param = "p_dhat", de = n, dp = dim_est,
            err = Matrix::norm(abar - p_dhat, "f")),
        list(est = "abar", param = "p_dhat_p", de = n, dp = d_p,
            err = Matrix::norm(abar - p_dhat_p, "f")),
        list(est = "phat", param = "p", de = dim_est, dp = n,
            err = Matrix::norm(phat - p, "f")),
        list(est = "phat", param = "p_dhat", de = dim_est, dp = dim_est,
            err = Matrix::norm(phat - p_dhat, "f")),
        list(est = "phat", param = "p_dhat_p", de = dim_est, dp = d_p,
            err = Matrix::norm(phat - p_dhat_p, "f")))

}

# norm(compute_phat(alist, dim = dim_est)$phat -
#     compute_phat(alist, dim = dim_est-1)$phat, "f")

# norm(low_rank_approx(abar, dim_est, FALSE) - 
#     low_rank_approx(p, dim_est, FALSE), "f")


# 100 replicates
nmc <- 100
# for each m in 1,5,10, 20, 50
mrange <- cross_df(list(mc = 1:nmc, m = c(1, 2, 5, 10, 20, 50, 100)))

err_sample_df <-  mrange %>% group_by(mc) %>%
    mutate(res = map(m, function(m){
        print(mc[1])
        # sample m graphs
        sample_df <- data_df %>% ungroup() %>% sample_n(m)

        # compute all the errors
        err_df <- try(sample_df %>%
            select(adj) %>% unlist %>%
            compare_to_latent(P, p_dec, p_dhat_p, d_p))

        # get rid of the adjacencies and store everything together
        sample_df %>% select(-adj) %>%
            nest(-atlas, -size, .key = "sample") %>%
            mutate(err = list(err_df))
    }))


save(err_sample_df, file = paste0("../../Data/test_latent_",
    aoi, "_", doi, ".RData"))

# Get the stuff we want
err_sample_df <- err_sample_df %>% unnest() %>% unnest(err)

g <- err_sample_df %>% filter(m < 400) %>%
    mutate(param = recode(param,
        p = "P", p_dhat_p = "P d*", p_dhat = "P dhat")) %>%
    mutate(err = err^2) %>%
        group_by(mc, m, param) %>%
    mutate(`relative efficiency` = ifelse(param == "P",
        err/err[est == "abar"],
        err/err[est == "abar"])) %>%
    filter(est == "phat") %>%
    ggplot(aes(x = m, y = `relative efficiency`))+
    stat_summary(fun.data = "mean_cl_boot") +
    stat_summary(fun.y = "mean", geom = "line") +
    facet_wrap(. ~ param) + geom_hline(yintercept = 1, linetype = 2)
g
ggsave(paste0("../../ieee_tmi_response/estimating_latent_",
    aoi, "_", doi, ".pdf"),
    width = 7, height = 3)

g + scale_x_log10() + scale_y_log10() +
    annotation_logticks()
ggsave(paste0("../../ieee_tmi_response/estimating_latent_",
    aoi, "_", doi, "_log.pdf"),
    width = 7, height = 3)

label_fun <- function(x){
    data.frame(label = format(mean(10^x), digits = 2), y = mean(x) - .5)
}

err_sample_df %>%
    mutate(err = err^2) %>%
    filter( (param == "p" & est == "abar") |
        (param == "p_dhat_p" & est == "phat")) %>%
    mutate(which = recode(param, p = "Abar estimating P",
        p_dhat_p = "Phat estiating P d*")) %>%
    ggplot(aes(x = m, y = err, color = which, group = which)) +
    stat_summary(fun.data = "mean_cl_boot", shape = 5) +
    stat_summary(fun.y = "mean", geom = "line") +
    stat_summary(fun.data = label_fun, geom = "text", color = "black") +
    scale_color_discrete(" ") + scale_x_log10() + scale_y_log10() +
    theme(legend.position="bottom") +
    annotation_logticks()

ggsave(paste0("../../ieee_tmi_response/estimating_latent_err_compare_",
    aoi, "_", doi, ".pdf"),
    width = 5, height = 5)

err_sample_df %>%
    mutate(err = err^2) %>%
    filter( (param == "p" & est == "abar") |
        (param == "p_dhat_p" & est == "phat")) %>%
    group_by(mc, m) %>% mutate(re = err/err[est == "abar"]) %>%
    filter(est == "phat") %>% group_by(m) %>% summarise(mre = mean(re))


