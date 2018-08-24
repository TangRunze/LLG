#! /usr/bin/env Rscript

#SBATCH -o /n/home15/dsussman/log/out-%a.txt
#SBATCH -e /n/home15/dsussman/log/err-%a.txt
#SBATCH -p stats
#SBATCH --mem-per-cpu=1500
#SBATCH -t 300
#SBATCH -a 1-200

# Set total jobs
total_jobs <- 200

library(methods)

# Get job ID if its available
if ( Sys.getenv("SLURM_JOB_ID") != "" ){
    job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    run_id <- Sys.getenv("SLURM_JOB_ID")
    data_dir <- "/n/regal/airoldi_lab/sussman/neurodata/"
    total_jobs <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
    server <- TRUE
} else {
    job_id <- 99
    run_id <- "local"
    server <- FALSE
    data_dir <- "/Volumes/Other/Data/neurodata_fmri/"
    total_jobs <- 1
}

print(run_id)
print(job_id)
print(total_jobs)


source("main_functions.R")
source("parse_neurodata.R")

save_dir <- paste0(data_dir, "results/")

aoi <- "DS01216_res-2x2x2"
aoi <- "CPAC200_res-2x2x2"
doi <- "SWU4"

# Load P and its decompositions so we don't have to do that everytime
load(paste0(data_dir, doi, "_", aoi, "_P.RData"))

# Load data
data_df <- df_fm(atlas_fm_df, fmri_scan_list) %>%
    filter(atlas == aoi, dataset == doi) %>%
    append_url_fn() %>%
    mutate(idx = seq(n())) %>%
    group_by_all() %>%
    select(-url)


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
nmc <- 5
# for each m in 1,5,10, 20, 50
mrange <- cross_df(list(mc = 1:nmc, m = c(1, 2, 5, 10, 20, 50, 100)))

set.seed(1000 + job_id)

read_edgelist <- function(fn){
    read_csv(fn, col_names = c("i", "j", "weight")) %>%
        from_data_frame(directed = FALSE)
}

err_sample_df <-  mrange %>% group_by(mc) %>%
    mutate(res = map(m, function(m){
        print(mc[1])
        # sample m graphs
        sample_df <- data_df %>% ungroup() %>%
            sample_n(m) %>%
            mutate(adj = list(read_graph(
                paste0(data_dir, "edgelist/", fn, ".edgelist"),
                format = "edgelist")[])) %>%
            ungroup() %>%
            select(-url, -dataset)

        # compute all the errors
        err_df <- try(sample_df %>%
            select(adj) %>% unlist %>%
            compare_to_latent(P, p_dec, p_dhat_p, d_p))

        # get rid of the adjacencies and store everything together
        sample_df %>% select(-adj) %>%
            nest(-atlas, -size, .key = "sample") %>%
            mutate(err = list(err_df))
    }))

err_sample_df$job_id <- job_id
err_sample_df$run_id <- run_id


save(err_sample_df,
    file = paste0(save_dir,
        paste(aoi, doi, run_id, job_id,
            "_test_latent_fmri.RData", sep = "_")))

pattern <- paste("test_latent", aoi, doi, run_id,
    "*fmri.RData", sep = "_")
saved <- list.files(save_dir, pattern)

if ( length(saved) == total_jobs ){
    cat("\n")
    all_df <- saved %>%
        map_df(function(x){
            cat(x, "\r")
            load(x);
            err_sample_df
    })
    cat("\n")
    save(all_df, file=paste0(wd,"/",pattern, "_ALL.RData"))
}
