library(tidyverse)
library(stringr)
library(igraph)

source("../../data/neurodata_list.R")

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

read_neurodata_graph <- function(d){
    read_graph(output_url(d), format="graphml")
}

load_all_graphs <- function(d){
    d %>% by_row(read_neurodata_graph,.to="igraph")
}

load_all_atlas <- function(d,da_df){
        d %>% filter(dataset==da_df$dataset) %>% 
        mutate(atlas=da_df$atlas[1]) %>%
        by_row(read_neurodata_graph,.to="igraph")   
}

igraph_to_df <- function(g){
    g <- as.matrix(g[])
    n <- nrow(g)
    expand.grid(i=1:n, j=1:n) %>% mutate(g=c(g)) %>% filter(i>j)
}

df_to_mat <- function(d,n){
    m <- matrix(0,n,n)
    m[d$i,d$j] <- d$P
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

neurodata_df <- neurodata_df %>% filter(atlas %in% atlas_keep)

all_desikan <- neurodata_df %>% filter(atlas=="desikan") %>% load_all_graphs()
save(all_desikan,file="/Volumes/Other/Data/neurodata_dtmri/all_desikan.RData")


a <- load_all_graphs(neurodata_df[1,],dataset_atlas_df %>% 
        filter(dataset=="HNU1"))

b <- all_desikan %>% group_by(subject,scan) %>%
    mutate(g=list(igraph_to_df(igraph[[1]]))) %>%
    select(-igraph) %>%
    unnest()