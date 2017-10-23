library(tidyverse)
library(ggpubr)
source("function_collection.R")
source("getElbows.R")

load_data <- FALSE
if( load_data ){

    data_dir <- "../../Data/mouse/"

    # Load ROI adjacency
    spat_adj <- (R.matlab::readMat(paste0(data_dir,
        "neighbors.mat"))$myneighbors>0) %>% which(arr.ind = TRUE) %>% 
            as_tibble() %>% 
            filter(row > col)

    # Load Connectivity
    conn_ll <- R.matlab::readMat(paste0(data_dir,
        "CHASSCONN148/connLL148.mat"))$connLL
    conn_lr <- R.matlab::readMat(paste0(data_dir,
        "CHASSCONN148/connLR148.mat"))$connLR
    conn_rl <- R.matlab::readMat(paste0(data_dir,
        "CHASSCONN148/connRL148.mat"))$connRL
    conn_rr <- R.matlab::readMat(paste0(data_dir,
        "CHASSCONN148/connRR148.mat"))$connRR

    # Bind into one 296x296 matrix and take the log
    conn_adj <- log(rbind(cbind(conn_ll,conn_lr),cbind(conn_rl,conn_rr))+1)

    # work with L1 for now
    label_name <- "L1"
    # read roi info
    roi_label_df <- readxl::read_xlsx(paste0(data_dir,
        "CHASSCONN148/labels_vals_plus2.xlsx")) %>% 
        filter(L1 != "exclude") %>% 
        select_(superstructure=label_name) %>%
        crossing(hemisphere=c("left","right"),.) %>%
        mutate(label = paste(hemisphere,superstructure)) 
    roi_label_ch <- roi_label_df$label
    roi_label <- roi_label_ch %>% as.factor() %>% as.numeric()


    #' Get xHat using diagonal augmentation and ZG for dim select
    get_x_hat <- function(m){

        n <- nrow(m)
        isSVD <- 1

        m_diag_aug <- diag_aug(m)
        nElbow <- 3
        evalVec <- ase(m_diag_aug, ceiling(n*3/5))[[1]]
        dZG <- getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
        A.ase <- ase(m_diag_aug, dZG, isSVD)
        x_hat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))
        x_hat[, 1:dZG]
    }

    x_hat <- get_x_hat(conn_adj)

    # compute data frame with all pairwise distances
    all_pair_diff <- dist(x_hat) %>%
        broom::tidy() %>%
        as_tibble %>%
        rename(i = item1, j = item2)

}

#' a function to randomly switch the labels for 2 pairs of ROIs
#' which cross the same boundary between labels
permute_roi <- function(spat_adj,roi_label){
    # keep all adjacent pairs with different labels
    spat_adj_diff <- spat_adj %>%
        mutate(lr = roi_label[row], lc = roi_label[col]) %>%
        filter(lr != lc)

    while(TRUE){
        # sample first pair
        s1 <- sample_n(spat_adj_diff, 1)
        # sample a different pair with the same labels
        s2_options <- spat_adj_diff %>% 
            filter(lr == s1$lr, lc== s1$lc, row != s1$row, col != s1$col)
        if(nrow(s2_options) == 0){
            next # we found an impossible situation
        }

        s2 <- spat_adj_diff %>% 
            filter(lr == s1$lr, lc== s1$lc, row != s1$row, col != s1$col) %>%
            sample_n(1)
        break
    }

    new_roi_label <- roi_label
    new_roi_label[c(s1$row,s2$row)] <- s1$lc
    new_roi_label[c(s1$col,s2$col)] <- s1$lr

    new_roi_label
}

# repeatedly permute the labels n_permute times
permute_roi_multi <- function(spat_adj, roi_label, n_permute){
    if(n_permute == 0){
        return(roi_label)
    }
    for(i in 1:n_permute){
        roi_label <- permute_roi(spat_adj, roi_label)
    }
    roi_label
}

# compute the difference between the mean distance within and the 
# mean difference between labels
cross_label_distance_metric <- function(all_pair_diff, roi_label) {
    all_pair_diff %>% ungroup() %>%
        mutate(same_label = roi_label[i] == roi_label[j]) %>%
        group_by(same_label) %>% 
        summarize(mean_diff = mean(distance)) %>%
        summarize(dmd = mean_diff[same_label]-mean_diff[!same_label]) %>%
        as.numeric()
}


compute <- TRUE
if( compute ){

# the original metric
metric_0 <- cross_label_distance_metric(all_pair_diff,
    roi_label)
system.time({

nmc <- 1000
n_perm_max <- 10

# permute and compute cross label distances
res_df <- crossing(mc = 1:nmc, n_perm = 1:n_perm_max) %>% 
    group_by(mc, n_perm) %>% 
    mutate(metric = 
        cross_label_distance_metric(all_pair_diff,
            permute_roi_multi(
                spat_adj, roi_label, n_perm)))
})
}else{
    load("temp.RData")
}




#############################################
############### Violin Plots ################
#############################################

# Summary stats with p-values
p_val <- res_df %>% group_by(n_perm) %>% 
    summarize(p_val = mean(metric<metric_0))


annotate_y <- with(res_df,min(metric)-(max(metric)-min(metric))*.05)


gg_violinT <- 
    ggplot(data = res_df, 
        aes(x=factor(n_perm), y=metric))+
    geom_violin(draw_quantiles = T, show.legend = FALSE)+
    geom_hline(yintercept = metric_0, linetype = 2) +
    scale_linetype_manual(name = "true lobe assignment",
        values = "dashed", labels = "") +
    guides(fill=FALSE)+
    theme(legend.position="bottom")+
    annotate("text", x = 1:n_perm_max, y = annotate_y, 
    label = paste0("p=", p_val$p_val))+
    labs(title = "", x = "number of flips",
        y = "T(X, l)", fill = "")

#############################################



# Make specrtal embedding data frame

ase_df <- data.frame(id=1:296,roi_label_df,x_hat) %>% 
    as_tibble %>%
    mutate(roi=as.numeric(as.factor(label)))






#############################################
############ X2 v. X4 w/QDA #################
#############################################

# Compute QDA classifier
qda8 <- ase_df %>%
    MASS::qda(label~X2+X4,data=.)


# Make a grid and classify all points in the grid
qda_class_df <- crossing(
    X2=seq(-3.3,4,.01),
    X4=seq(-3,2.7,.01))
qda_pred <- predict(qda8,qda_class_df)
qda_class_df$label  <- predict(qda8,qda_class_df)$class

# Add in a ratio between biggest and second biggest
# posteriors to help determin boundaries
qda_class_df$post_ratio <- qda_pred$post %>%
    apply(1,function(x){
        r <- sort(x,decreasing=T); 
        r[2]/r[1]
    })

# Separate hemisphere and superstructure
qda_class_df <- qda_class_df %>%
    separate(label,into=c("hemisphere","superstructure"),
            sep=" ",extra="merge",remove=FALSE)

qda_mean_df <- tibble(label=row.names(qda8$means)) %>% 
    bind_cols(as_tibble(qda8$means))

post_ratio_threshold <- 0.925
Lhat <- round(1 - mean(predict(qda8,ase_df)$class == ase_df$label),2)

gg_class <- ase_df %>% group_by(id) %>% 
    ggplot(aes(x = X2, y = X4,
        group=label, color = superstructure, 
        fill=superstructure, alpha=hemisphere,
        linetype=hemisphere, shape=hemisphere)) +
    geom_raster(data=qda_class_df) +
    geom_raster(data=qda_class_df %>% 
            filter(post_ratio>post_ratio_threshold),fill="black",alpha=0.5)+
    geom_point(alpha=.6) +
    stat_ellipse(aes(fill=superstructure),
        level = 0, alpha=1,geom="point",
        size=3,color="black")+
    scale_shape_manual(values = c(21,22))+
    scale_alpha_discrete(range=c(0.4,0.3))
    #  +
    # annotate("text",x=1,y=-2.8,
    #     label=paste0("training error=",Lhat),size=4)
#############################################



#############################################
############ X1--X4 point plot ##############
#############################################
gg_4dim <- ase_df %>% 
    select(-X5, -X6, -X7, -roi) %>%
    gather(dimension,value,
        -label, -hemisphere, -superstructure, -id) %>%
    ggplot(aes(x = id,y = value, group = label,
        color = superstructure, shape = hemisphere)) +
    geom_point(alpha=.3) + 
    geom_smooth(method="lm", formula = y ~ 1) +
    scale_shape_manual(values = c(21,22))+
    facet_wrap(~dimension,scales="free_y")
#############################################



#############################################
################### X2 plot #################
#############################################
ase_df %>%
    ggplot(aes(x=id,y=X2,group=label,color=label)) +
    geom_point(alpha=1) +
    geom_smooth(method="lm", formula = y ~ 1)
#############################################




ase_plot <- ggarrange(gg_4dim, gg_class,
    common.legend = TRUE,
    ncol = 1,
    nrow = 2,
    legend = "bottom")
ggexport(ase_plot,
    filename="../../figure/mouse_connectome_ase.pdf",
    height = 8,
    width = 6,
    pointsize = 10)


ggexport(gg_violinT,
    filename="../../figure/mouse_connectome_violin.pdf",
    height = 3,
    width = 6,
    pointsize = 10)