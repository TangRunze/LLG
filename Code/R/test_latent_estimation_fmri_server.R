#! /usr/bin/env Rscript

#SBATCH -o /n/home15/dsussman/log/out-%a.txt
#SBATCH -e /n/home15/dsussman/log/err-%a.txt
#SBATCH -p stats
#SBATCH --mem-per-cpu=3200
#SBATCH -t 300
#SBATCH -a 1-2

# Set total jobs
total_jobs <- 200

library(methods)

# Get job ID if its available
if ( Sys.getenv("SLURM_JOB_ID") != "" ){
    job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    run_id <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
    data_dir <- "/n/regal/airoldi_lab/sussman/neurodata/"
    total_jobs <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
    server <- TRUE
} else {
    job_id <- 99
    run_id <- "local"
    server <- FALSE
    data_dir <- "/n/regal/airoldi_lab/sussman/neurodata/"
    total_jobs <- 1
}

print(run_id)
print(job_id)
print(total_jobs)


source("main_functions.R")
source("parse_neurodata.R")

save_dir <- paste0(data_dir, "results/")

aoi <- "DS01216_res-2x2x2"
# aoi <- "CPAC200_res-2x2x2"
doi <- "SWU4"

# Load P and its decompositions so we don't have to do that everytime
load(paste0(data_dir, doi, "_", aoi, "_P.RData"))

# Load data
data_df <- df_fm(atlas_fm_df, fmri_scan_list) %>%
    filter(atlas == aoi, dataset == doi) %>%
    append_url_fn() %>%
    mutate(idx = seq(n())) %>%
    group_by_all()


# compute P = average of all
m_all <- nrow(data_df)
n <- nrow(P)


compare_to_latent <- function(alist, p, p_dec, p_dhat_p, d_p){
    m <- length(alist)
    abar <- Reduce("+", alist) / m

    # Compute estimate
    phat <- compute_phat(alist)
    dim_est <- phat$dim
    phat  <- phat$phat

    p_dhat <- compute_phat(alist, abar = p, dim = dim_est)$phat

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
nmc <- 1
# for each m in 1,5,10, 20, 50
mrange <- cross_df(list(mc = 1:nmc, m = c(1, 2, 5, 10, 20, 50)))

set.seed(1000 + job_id)

read_edgelist <- function(fn){
    read_delim(fn, col_names = c("i", "j", "weight"),
            delim = " ", progress = FALSE, col_types = "ddd") %>%
        graph_from_data_frame(directed = FALSE)
}

err_sample_df <-  mrange %>% group_by(mc) %>%
    mutate(res = map(m, function(m){
        print(c(mc[1], m))
        gc()
        # sample m graphs
        sample_df <- data_df %>% ungroup() %>%
            sample_n(m) %>%
            group_by_all() %>%
            mutate(adj = list(read_edgelist(
                paste0(data_dir, "edgelist/",
                    fn, ".edgelist"))[])) %>%
            ungroup() %>%
            select(-url)

        # compute all the errors
        err_df <- sample_df %>%
            select(adj) %>% unlist %>%
            compare_to_latent(P, p_dec, p_dhat_p, d_p)

        # get rid of the adjacencies and store everything together
        sample_df %>% select(-adj) %>%
            nest(-atlas, -dataset, .key = "sample") %>%
            mutate(err = list(err_df))
    }))

err_sample_df$job_id <- job_id
err_sample_df$run_id <- run_id


save(err_sample_df,
    file = paste0(save_dir,
        paste(aoi, doi, run_id, job_id,
            "test_latent_fmri.RData", sep = "_")))

pattern <- paste(aoi, doi,
    "*_test_latent_fmri.RData", sep = "_") %>% glob2rx()
saved <- list.files(save_dir, pattern)

if ( length(saved) == total_jobs ){
    cat("\n")
    all_df <- saved %>%
        map_df(function(x){
            cat(x, "\r")
            load(paste0(save_dir, x))
            err_sample_df
    })
    cat("\n")
    save(all_df, file=paste0(save_dir,
        paste(aoi, doi, "all",
            "test_latent_fmri.RData", sep = "_")))
}

process <- function(){

fn <- "/Volumes/Other/Data/neurodata_fmri/test_latent_SWU4_CPAC200_res-2x2x2_51637946_fmri_ALL.RData"
fn <- "~/Dropbox (Personal)/Manuscript/LLG/Data/DS01216_res-2x2x2_SWU4_all_test_latent_fmri.RData"
load(fn)
unnest(all_df, res)$err[1]
err_sample_df <- all_df %>% unnest(res) %>% unnest(err)

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

g + scale_x_log10() + scale_y_log10() +
    annotation_logticks()





err_sample_df %>%
    mutate(err = err^2) %>%
    filter( (param == "p" & est == "abar") |
        (param == "p_dhat_p" & est == "phat")) %>%
    mutate(which = recode(param, p = "Abar estimating P",
        p_dhat_p = "Phat estiating P d*")) %>%
    ggplot(aes(x = m, y = err, color = which, group = which)) +
    stat_summary(fun.data = "mean_cl_boot", shape = 5) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_color_discrete(" ") + scale_x_log10() + scale_y_log10() +
    theme(legend.position="bottom") +
    annotation_logticks()
}
