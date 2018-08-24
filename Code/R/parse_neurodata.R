library(tidyverse)
library(stringr)
library(igraph)

source("neurodata_list.R")

parse_file_name <- function(s){
    d <- str_split(s,"_")[[1]]
    data.frame(dataset=d[1],subject=d[2],scan=d[3]) %>%
        as_tibble()
}

output_url <- function(d){
    fn <- paste0(with(d,
            paste(dataset,subject,scan,"DTI",atlas,sep="_")),
            ".graphml")
    with(d, paste("http://openconnecto.me/mrdata/share/dti/ndmg_v0011",
        dataset,atlas,fn,sep="/"))
}

output_url_fm <- function(d){

    fmri_scan_list
}


read_neurodata_graph <- function(d){
    read_graph(output_url(d), format = "graphml")
}

load_all_graphs <- function(d){
    d %>% by_row(read_neurodata_graph,.to="igraph")
}

load_all_atlas <- function(d, da_df){
        d %>% filter(dataset==da_df$dataset) %>% 
        mutate(atlas=da_df$atlas[1]) %>%
        by_row(read_neurodata_graph,.to="igraph")   
}

igraph_to_df <- function(g){
    g <- as.matrix(g[])
    n <- nrow(g)
    expand.grid(i = 1:n, j = 1:n) %>% mutate(g = c(g)) %>% filter(i > j)
}

df_to_mat <- function(d,n){
    m <- matrix(0,n,n)
    m[d$i,d$j] <- d$P
}

load_graphml <- function(dataset, subject, scan, atlas, size, data_dir){
    fn <- paste0(data_dir, dataset,
        "sub-", subject,
        "_ses-", scan,
        "_dwi_", atlas, ".graphml")
    tibble(adj = list(try(igraph::read_graph(fn, format = "graphml")[])))
}

neurodata_df <- neurodata_scan_list %>%
    map_df(parse_file_name) %>% 
    merge(atlas_df) %>% 
    as_tibble()


atlas_keep <- c(
    "AAL",
    "CPAC200",
    "HarvardOxford",
    "JHU",
    "Talairach",
    "desikan")

neurodata_df <- neurodata_df %>% filter(atlas %in% atlas_keep) %>%
    mutate(atlas = as.character(atlas))

df_fm <- function(atlas_list, scan_list){
    tibble(
        dataset = rep(names(scan_list),
            scan_list %>% map(length)),
        scan = unlist(fmri_scan_list)) %>%
        separate(scan, c("subject", "session"), sep = "_") %>%
        tidyr::crossing(atlas_list)

}

url_fm <- function(d){
    # Here is an example url
    # http://mrneurodata.s3.amazonaws.com/data/fmri/
    # SWU4/ndmg_0-0-1f/func/connectomes/CPAC200_res-2x2x2/
    # sub-0025629_ses-1_bold_CPAC200_res-2x2x2_measure-correlation.gpickle
    base_url <- "http://mrneurodata.s3.amazonaws.com/data/fmri/"
    url <- paste0(base_url, 
        d$dataset, "/ndmg_0-0-1f/func/connectomes/",
        d$atlas, 
        "/sub-", d$subject, "_ses-", d$session, "_bold_",
        d$atlas, "_measure-correlation.gpickle")
}

fn_fm <- function(d){
    paste0(paste(
        d$dataset, d$atlas, d$subject, d$session, sep = "_"))
}

append_url_fn <- function(df){
    df %>% mutate(url = url_fm(.), fn = fn_fm(.))
}

download_fm <- function(url, fn, data_dir, ...){
    print(fn)
    download.file(url,
        paste0(data_dir, fn, ".gpickle"), method = "auto")
}

# data_dir <- "/Volumes/Other/Data/neurodata_fmri/"
# # data_dir <- "~/Dropbox/Data/neurodata_fmri/"
# neuro_fm_df <- df_fm(atlas_fm_df, fmri_scan_list) %>%
#     filter(atlas == "DS01216_res-2x2x2", dataset == "SWU4") %>%
#     append_url_fn() %>%
#     pwalk(download_fm, data_dir = data_dir)

# dir("/Volumes/Other/Data/neurodata_fmri/graphml/")
# g <- read_graph("/Volumes/Other/Data/neurodata_fmri/graphml/SWU4_DS01216_res-2x2x2_0025629_2.graphml", format = "graphml")

# neuro_fm_df <- neuro_fm_df %>%
#     group_by_all() %>%
#     mutate(g = list(read_graph(
#         paste0(data_dir, "graphml/", fn, ".graphml"),
#         format = "graphml")))


# all_desikan <- neurodata_df %>% filter(atlas=="desikan") %>% load_all_graphs()
# save(all_desikan,file="/Volumes/Other/Data/neurodata_dtmri/all_desikan.RData")


# a <- load_all_graphs(neurodata_df[1,],dataset_atlas_df %>% 
#         filter(dataset=="HNU1"))

# b <- all_desikan %>% group_by(subject,scan) %>%
#     mutate(g=list(igraph_to_df(igraph[[1]]))) %>%
#     select(-igraph) %>%
#     unnest()