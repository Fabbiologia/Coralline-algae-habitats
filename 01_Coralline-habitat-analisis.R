## Script related to the paper ##
## "Coralline algae habitats: ##
## an environmental and microstructural ##
## assessment to understand microbial niches" ##
## Creator Fabio Favoretto: favoretto.fabio@gmail.com ##


# Loading libraries -------------------------------------------------------


library(RCurl)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(vegan)
library(ggthemes)
library(mvpart)
library(MVPARTwrap)




# Loading data ------------------------------------------------------------


level_5 <- read.csv(url("https://raw.githubusercontent.com/Fabbiologia/Coralline-algae-habitats/master/data/level_5.csv"))
halophilas <- read.csv(url("https://raw.githubusercontent.com/Fabbiologia/Coralline-algae-habitats/master/data/halophilas.csv"))
env <- read.csv(url("https://raw.githubusercontent.com/Fabbiologia/Coralline-algae-habitats/master/data/env.csv"))





# Heatmap Taxa level 5 ----------------------------------------------------


# Matrix creation 

mtx2 <- as.matrix(level_5[2:ncol(level_5)]) # this transform the data frame to a matrix (needed for the heatmap)


rownames(mtx2) <- level_5$index # this apply rownames to the matrix (same as sample names)
tmtx <- t(mtx2)
tmtx <- tmtx[,1:12]

# heatmap plotting

col_fun <- colorRamp2(c(-3, 0, 3), c("dodgerblue4", "white", "firebrick")) #select colour ramp

Heatmap(t(scale(t(tmtx))), 
        name = "Scale", 
        cluster_columns = FALSE,
        col = col_fun,
        clustering_distance_rows = "pearson",
        column_split = factor(c(rep("Winter",6),rep("Summer", 6)), levels = c("Winter", "Summer")),
        column_title_gp = gpar(fill = c("#D0EAE1", "#EBB6AC")),
        row_dend_reorder = TRUE,
        row_names_gp = gpar(fontsize = 5), 
        row_km = 2) 




# NDMS model --------------------------------------------------------------

alg <- level_5[1:12,] # selecting only algae samples

# wrangling data:
spe <- alg %>% 
        select(-c(index))

rownames(spe) <- alg$index

grouping_info <- data.frame(
        sample = rownames(spe),
        season = c(rep("Winter", 6), rep("Summer", 6)),
        site = c(rep("Pool_1", 3), rep("Pool_2", 3), rep("Pool_1", 3), rep("Pool_2", 3)))


#Get MDS stats

sol <- metaMDS(spe,distance = "bray", trymax = 50)

plot.new()
ord <- ordiellipse(sol, as.factor(grouping_info[,2]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

NMDS <- data.frame(x=sol$point[,1],y=sol$point[,2],Season=as.factor(grouping_info[,2]),Site=as.factor(grouping_info[,3]))

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
        theta <- (0:npoints) * 2 * pi/npoints
        Circle <- cbind(cos(theta), sin(theta))
        t(center + scale * t(Circle %*% chol(cov)))
}

# Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Season)){
        if(g!="" && (g %in% names(ord))){
                
                df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Season==g,],
                                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                              ,Season=g))
        }
}
# Generate mean values from NMDS plot grouped on Seasons
NMDS.mean <- aggregate(NMDS[,1:2],list(group=NMDS$Season),mean)

shape_values<-seq(1,11) #graphical parameters


# Standardizing environmental variables
env.z <- as.data.frame(decostand(env, method="standardize"))

# Fitting environmental variables
env.fit <- envfit(sol, env.z, permutations = 9999)

(env.fit$vectors)

env_sel <- env.z %>% 
        select(temperature, N_NH4, N_NO3, N_NO2, Si_SiO2, DIN, TN, TNtoTP)
#Get site information
df<-scores(sol,display=c("sites"))

#Add grouping information
df<-data.frame(df,Type=grouping_info[rownames(df),1])


vec.sp<-envfit(sol$points, env_sel, perm=9999)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)



ggplot(data=NMDS,aes(x,y,colour=Season))+ 
        annotate("text",x=NMDS.mean$x,y=NMDS.mean$y+0.4,label=NMDS.mean$group,size=4, fontface = 2)+ 
        geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=.1, linetype=1)+
        geom_point(aes(shape=Site))+
        labs(caption = "Stress = 0.02")+
        scale_color_brewer(palette = "Set2", direction = -1)+
        geom_segment(data=vec.sp.df,
                     aes(x=0,xend=MDS1,y=0,yend=MDS2), 
                     arrow = arrow(length = unit(0.1, "cm")),
                     colour="grey90") + 
        ggrepel::geom_text_repel(data = vec.sp.df, 
                                 aes(x=MDS1, y=MDS2, label=species), 
                                 size=4, inherit.aes = FALSE)+
        theme_classic()+
        theme(axis.line = element_blank(), 
              axis.ticks = element_blank(),
              axis.text = element_blank(), 
              axis.title = element_blank(), 
              panel.border = element_rect(fill = NA), 
              legend.position = "")

tiff(filename = "NMDS.tif", width = 8, height = 5, units = "in", res = 300)

ggplot(data=NMDS,aes(x,y,colour=Season))+ 
        annotate("text",x=NMDS.mean$x,y=NMDS.mean$y+0.4,label=NMDS.mean$group,size=4, fontface = 2)+ 
        geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=.1, linetype=1)+
        geom_point(aes(shape=Site))+
        labs(caption = "Stress = 0.02")+
        scale_color_brewer(palette = "Set2", direction = -1)+
        geom_segment(data=vec.sp.df,
                     aes(x=0,xend=MDS1,y=0,yend=MDS2), 
                     arrow = arrow(length = unit(0.1, "cm")),
                     colour="grey90") + 
        ggrepel::geom_text_repel(data = vec.sp.df, 
                                 aes(x=MDS1, y=MDS2, label=species), 
                                 size=4, inherit.aes = FALSE)+
        theme_classic()+
        theme(axis.line = element_blank(), 
              axis.ticks = element_blank(),
              axis.text = element_blank(), 
              axis.title = element_blank(), 
              panel.border = element_rect(fill = NA), 
              legend.position = "")
dev.off()

# Multivariate Regression Tree --------------------------------------------

# species transformation 
spe.hel <- decostand(spe, method="hellinger") 

# Create the regression tree
env_sel <- env %>% 
        select(temperature, N_NO3, N_NO2, TN, Si_SiO2, DIN, Chl_a)

var <- env[1:13]
spe.hel <- decostand(spe, method="hellinger") # you can also use method="hell" 

spe.hel[is.na(spe.hel)] <- 0


# Interactive group selection WARNING!! need to clic esc to exit the interactive visualization

doubs.mrt <- mvpart(as.matrix(spe.hel) ~. ,env_sel,
                    legend=FALSE, margin=0.01, cp=0, xv="pick",
                    xval=nrow(spe), xvmult=100, which=4,
                    pca = TRUE)


rpart.pca(doubs.mrt, interact=TRUE)
summary(doubs.mrt)

# Find discriminant species with MRT results
doubs.mrt.wrap<-MRT(doubs.mrt,percent=10,species=colnames(spe.hel))
summary.MRT(doubs.mrt.wrap)

# Extract indval p-values
doubs.mrt.indval<-indval(spe.hel,doubs.mrt$where)
doubs.mrt.indval$pval

# Extract indicator species of each node, with its indval
(doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)])
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]






# Heatmap halophilas ------------------------------------------------------

## data wrangling of halophilas dataset
mtx <- as.matrix(halophilas[2:ncol(halophilas)]) # this transform the data frame to a matrix (needed for the heatmap)
rownames(mtx) <- halophilas$index # this apply rownames to the matrix (same as sample names)
mtx <- mtx[1:12,]

selMtx<- vegan::decostand(mtx, method = 'pa') %>% 
        as.data.frame(.) %>% 
        colSums(.) %>% 
        as.data.frame(.) %>% 
        dplyr::filter(.>=2)

mtx2 <- as.matrix(halophilas[2:ncol(halophilas)]) # this transform the data frame to a matrix (needed for the heatmap)


rownames(mtx2) <- halophilas$index # this apply rownames to the matrix (same as sample names)

tmtx <- t(mtx)
tmtx <- as.data.frame((tmtx))
fullnames <- rownames(as.data.frame(tmtx))

tmtx$id <- fullnames
tmtx <- tmtx %>% 
        dplyr::filter(id %in% rownames(selMtx)) %>% 
        select(-id)

tmtx2 <- tmtx[,1:12]
ind <- apply(tmtx2, 1, function(x) all(x==0))
tmtx2 <- tmtx2[ !ind, ]

## end of data wrangling

# heatmap plotting 
col_fun <- colorRamp2(c(-3, 0, 3), c("dodgerblue4", "white", "firebrick")) # color ramp

Heatmap(t(scale(t(tmtx2))), 
        name = "Scale", 
        cluster_columns = FALSE,
        col = col_fun,
        clustering_distance_rows = "pearson",
        row_dend_reorder = TRUE,
        row_names_gp = gpar(fontsize = 8)) 




## END OF SCRIPT ==========

