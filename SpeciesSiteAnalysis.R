library(ggplot2)
library(stats) #glms
library(ggthemes)
library(car)
library(scales)
library(AER) #has dispersiontest() to check for overdispersion in poisson
library(tidyr) #gather
library(ggplot2)
library(colortools)
library(MASS)

df <- read.csv('Kleptoparasites_FuncEcol_data.csv')
df$site.elevation <- as.factor(df$site.elevation)

### get klepto presence / absence from counts ####
df$klepto.presence <- rep(0, times=length(df$nest))
rows <- which(df$klepto_all == 0) # keep as is
rows <- which(df$klepto_all > 0)
df$klepto.presence[rows] <- 1

# get column with species name
df$species <- rep(0, times=length(df$nest))
rows <- which(df$species.site == 'Adom400')
df$species[rows] <- 'domingo'
rows <- which(df$species.site=='Aexi400')
df$species[rows] <- 'eximius'
rows <- which(df$species.site=='Aexi1150')
df$species[rows] <- 'eximius'
rows <- which(df$species.site=='Aele1150')
df$species[rows] <- 'elegans'
rows <- which(df$species.site=='Aele1600')
df$species[rows] <- 'elegans'
rows <- which(df$species.site=='Agua1600')
df$species[rows] <- 'guacamayos'

# reorder factors
df$site.elevation <- factor(df$site.elevation, levels=c("400", "1150", "1600"), labels=c("400", "1150", "1600"))
df$species <- factor(df$species, levels=c("eximius", "domingo", "elegans", "guacamayos"), labels=c("eximius", "domingo", "elegans", "guacamayos"))

# include only social species at two lower elevations
df2 <- df[-which(df$site.elevation == '1600'),]
df2$klepto_per_capita <- (df2$klepto_all) / (df2$Anelosimus_adults_subadults)

# df3 <- df2[-which(df2$klepto_all == 0),]

########### Figure 2 analyses ###############
# Figure 2a
z <- glm(klepto.presence ~ Log10tangleVol_cm3 + site.elevation + Log10tangleVol_cm3:site.elevation, family = binomial, data = df)
anova(z, test = "Chisq")
Anova(z, type = "II")

F2a <- ggplot(df, aes(x=tangleVol_cm3, y=klepto.presence, linetype = factor(site.elevation), colour = factor(site.elevation))) + 
        geom_point(size = 4) + 
        stat_smooth(method="glm", method.args = list(family="binomial"), se=F, size = 1.5) +
        scale_colour_manual(values = c("#000000", "#999999","#CCCCCC"))+
        scale_linetype_manual(values = c("solid", "dashed","twodash"))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "Kleptoparasite presence",
             colour = "Site Elevation",
             linetype = "Site Elevation")+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24),
              # theme(legend.position = 'none')
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width=unit(3,"line"))+
        theme(legend.justification=c(1.03,0), legend.position=c(1,0.16))
F2a

# Figure 2b
y <- glmer(klepto.presence ~ nest_elevation + (1|species), family = binomial, data = df)
summary(y)
Anova(y, type = "II")

F2b <- ggplot(df, aes(x=nest_elevation, y=klepto.presence)) + 
        geom_point(size=4) + 
        stat_smooth(method="glm", method.args = list(family="binomial"), se=F, size = 1.5, colour = 'black') +
        labs(x = "Nest elevation (m)", 
             y = "Kleptoparasite presence")+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))
F2b

multiplot(F2b, F2a)

########## Figure 3 analyses #############
#figure 3a analyses
p <- glm(klepto_all ~ site.elevation/species, data = df, family = poisson)
Anova(p, type = "II") 
dispersiontest(p,trafo=1, alternative = 'greater') #overdispersed

library(MASS)
p2 <- glm.nb(klepto_all ~ site.elevation/species, data = df)
Anova(p2, type = "II")

# Figure 3a ggplot
b = c(1,2, 5, 10, 25, 50)
F3a <- ggplot(df, aes(x=site.elevation, y=klepto_all, fill = species))+
        geom_boxplot()+
        theme_base()+
        scale_fill_manual(values = c("#00b979ff", "#554236","#d3ce3d", "#f1efa5"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#f7bb94ff", "#d3ce3d"))+
        scale_y_log10(breaks = b,
                    labels = b)+
        labs(x = "Site elevation",
             y = "# Kleptoparasites",
                fill = "Species", colour = "Species")+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24),
              legend.text=element_text(size=18),
              legend.title = element_text(size=22),
              legend.position=c(0.75,0.8))

F3a

#Fig3b analysis
p <- glm(KDPR_percent ~ site.elevation/species, data = df, family = gaussian)
Anova(p, type = "II")

F3b <- ggplot(df, aes(x=site.elevation, y=(KDPR_percent), fill = species))+
        geom_boxplot()+
        theme_base()+
        scale_fill_manual(values = c("#00b979ff", "#554236","#d3ce3d", "#f1efa5"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#f7bb94ff", "#d3ce3d"))+
        scale_x_discrete(name = "Site elevation") +
        scale_y_continuous(name = "% Performing Behaviours")+
        labs(fill = "Species", colour = "Species")+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24),
              legend.position='none')
F3b

multiplot(F3a, F3b, cols=2)

####### Figure 4 analyses ########
# Fig4a

df3$klepto_all %>% hist()
df3$tangleVol_cm3 %>%  hist()
df3$tangleVol_cm3 %>% log() %>% hist()

klepto.tangle <- glm(klepto_all ~ Log10tangleVol_cm3 + species.site + Log10tangleVol_cm3:species.site, 
                     data = df2, family = poisson(link='log'))
Anova(klepto.tangle, type = "II")
summary(klepto.tangle)

# pick apart Log10tangleVol_cm3:species.site
# Adom400
####currently with df3 (no zeroes)
rows <- which(df2$species.site=="Adom400")
Adom400 <- df2[rows,]
F2a.adom400 <- glm(klepto_all ~ Log10tangleVol_cm3, 
                   data = Adom400, family = poisson(link='log'))
Anova(F2a.adom400, type = "III")
summary(F2a.adom400)
p1 <- ggplot(Adom400, aes(x =(tangleVol_cm3), y = (klepto_all)))+
        geom_point() + 
        ggtitle("A. domingo 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=poisson(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites")+
        expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aexi400
rows <- which(df2$species.site == "Aexi400")
Aexi400 <- df2[rows,]
F2a.aexi400 <- glm(klepto_all ~ Log10tangleVol_cm3, 
                   data = Aexi400, family = poisson(link='log'))
Anova(F2a.aexi400, type = "III")
summary(F2a.aexi400)
p2 <- ggplot(Aexi400, aes(x =(tangleVol_cm3), y = (klepto_all)))+
        geom_point() + 
        ggtitle("A. eximius 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=poisson(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        expand_limits(y=c(0,70))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites")+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aexi1150
rows <- which(df2$species.site == "Aexi1150")
Aexi1150 <- df2[rows,]
F2a.aexi1150 <- glm(klepto_all ~ Log10tangleVol_cm3, 
                   data = Aexi1150, family = poisson(link='log'))
Anova(F2a.aexi1150, type = "III")
summary(F2a.aexi1150)
p3 <- ggplot(Aexi1150, aes(x =(tangleVol_cm3), y = (klepto_all)))+
        geom_point() + 
        ggtitle("A. eximius 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=poisson(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites")+
        expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aele1150
rows <- which(df2$species.site == "Aele1150")
Aele1150 <- df2[rows,]
F2a.aele1150 <- glm(klepto_all ~ Log10tangleVol_cm3, 
                    data = Aele1150, family = poisson(link='log'))
Anova(F2a.aele1150, type = "III")
p4 <- ggplot(Aele1150, aes(x =(tangleVol_cm3), y = (klepto_all)))+
        geom_point() + 
        ggtitle("A. elegans 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=poisson(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites")+
        expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

multiplot(p3, p2, p4, p1, cols=2)
## Figure 4a code ######
F4a <- ggplot(df2, aes(x =(tangleVol_cm3), y = (klepto_all), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size = 4) + 
        stat_smooth(method="glm", method.args = list(family=poisson(link='log')), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites",
             colour = "Species",
             linetype = "Site Elevation",
             shape = 'Species')+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))

F4a

##### Fig4b analyses #######
df3 <- df2[-which(df2$klepto_all == 0),]
kleptoDens.tangle <- glm(KleptoDens_dm3 ~ Log10tangleVol_cm3 + species.site + Log10tangleVol_cm3:species.site, 
                         data = df3, family = Gamma(link='log'))
Anova(kleptoDens.tangle, type = "II")

plot(kleptoDens.tangle)


rows <- which(df3$species.site=="Adom400")
Adom400 <- df3[rows,]
rows <- which(df3$species.site == "Aexi400")
Aexi400 <- df3[rows,]
rows <- which(df3$species.site == "Aexi1150")
Aexi1150 <- df3[rows,]
rows <- which(df3$species.site == "Aele1150")
Aele1150 <- df3[rows,]

# Adom400
F2b.adom400 <- glm(KleptoDens_dm3 ~ Log10tangleVol_cm3, 
                   data = Adom400, family = Gamma(link='log'))
Anova(F2b.adom400, type = "III")
summary(F2b.adom400)
p1 <- ggplot(Adom400, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3)))+
        geom_point() + 
        ggtitle("A. domingo 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}))+
        # expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))


# Aexi400
F2b.aexi400 <- glm(KleptoDens_dm3 ~ Log10tangleVol_cm3, 
                   data = Aexi400, family = Gamma(link='log'))
Anova(F2b.aexi400, type = "III")
summary(F2b.aexi400)
p2 <- ggplot(Aexi400, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3)))+
        geom_point() + 
        ggtitle("A. eximius 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}))+
        # expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aexi1150
F2b.aexi1150 <- glm(KleptoDens_dm3 ~ Log10tangleVol_cm3, 
                   data = Aexi1150, family = Gamma(link='log'))
Anova(F2b.aexi1150, type = "III")
summary(F2b.aexi1150)
p3 <- ggplot(Aexi1150, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3)))+
        geom_point() + 
        ggtitle("A. eximius 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}))+
        # expand_limits(y=c(0,70))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aele1150
F2b.aele1150 <- glm(KleptoDens_dm3 ~ Log10tangleVol_cm3, 
                    data = Aele1150, family = Gamma(link='log'))
Anova(F2b.aele1150, type = "III")
p4 <- ggplot(Aele1150, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3)))+
        geom_point() + 
        ggtitle("A. elegans 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}))+
        # expand_limits(y=c(0,50))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

multiplot(p3, p2, p4, p1, cols = 2)

######### Figure 4b code #######

df3$KleptoDens_dm3 %>% hist()
df3$tangleVol_cm3 %>% log() %>% hist()

F4b <- ggplot(df3, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size = 4) + 
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}),
             colour = "Species",
             linetype = "Site Elevation",
             shape = 'Species')+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))
F4b

###### Fig4c analysis #########
df3$klepto_per_capita %>% hist()
hist(1/df$Anelosimus_adults_subadults)

kleptoDens.per.cap <- glm(klepto_per_capita ~ log10(Anelosimus_adults_subadults) + species.site + log10(Anelosimus_adults_subadults):species.site, 
                          data = df3, family = Gamma(link='log'))
Anova(kleptoDens.per.cap, type = "II")
plot(kleptoDens.per.cap)


# Adom400
F2b.adom400 <- glm(klepto_per_capita ~ log10(Anelosimus_adults_subadults), 
                   Adom400, family = Gamma(link='log'))
Anova(F2b.adom400, type = "III")
summary(F2b.adom400)
p1 <- ggplot(Adom400, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita)))+
        geom_point() + 
        ggtitle("A. domingo 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita")+
        expand_limits(y=c(0,1))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))
p1
# Aexi400
F2b.aexi400 <- glm(klepto_per_capita ~ log10(Anelosimus_adults_subadults), 
                   data = Aexi400, family = Gamma(link='log'))
Anova(F2b.aexi400, type = "III")
summary(F2b.aexi400)
p2 <- ggplot(Aexi400, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita)))+
        geom_point() + 
        ggtitle("A. eximius 400") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita")+
        expand_limits(y=c(0,1))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aexi1150
F2b.aexi1150 <- glm(klepto_per_capita ~ log10(Anelosimus_adults_subadults), 
                    data = Aexi1150, family = Gamma(link='log'))
Anova(F2b.aexi1150, type = "III")
p3 <- ggplot(Aexi1150, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita)))+
        geom_point() + 
        ggtitle("A. eximius 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita")+
        expand_limits(y=c(0,1))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

# Aele1150
F2b.aele1150 <- glm(klepto_per_capita ~ log10(Anelosimus_adults_subadults), 
                    data = Aele1150, family = Gamma(link='log'))
Anova(F2b.aele1150, type = "III")
summary(F2b.aele1150)
p4 <- ggplot(Aele1150, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita)))+
        geom_point() + 
        ggtitle("A. elegans 1150") +
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=Gamma(link='log')), se=T, size = 1.5) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita")+
        expand_limits(y=c(0,1))+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))

multiplot(p3, p2, p4, p1, cols=2)


F4c <- ggplot(df3, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size = 4) + 
        stat_smooth(method="glm", method.args = list(family=Gamma(link='log')), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita",
             colour = "Species",
             linetype = "Site Elevation",
             shape = 'Species')+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))

F4c

# Fig4c
kleptoDens.per.cap <- glm(log10(klepto_per_capita) ~ log10(Anelosimus_adults_subadults) + species.site + log10(Anelosimus_adults_subadults):species.site, 
                          data = df3, family = gaussian)
Anova(kleptoDens.per.cap, type = "III")


# ######## Possible Fig 4d #######
F4d <- ggplot(df3, aes(x =(tangleVol_cm3), y = (klepto_all), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size=4) + 
        stat_smooth(method="glm", method.args = list(family=gaussian), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = "# Kleptoparasites",
             colour = "Species",
             linetype = "Site Elevation", 
             shape = "Species")+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))

F4d

## Fig4e
F4e <- ggplot(df3, aes(x =(tangleVol_cm3), y = (KleptoDens_dm3), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size=4) + 
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=gaussian(link='identity')), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = expression(Tangle~Volume~cm^{3}), 
             y = expression(Kleptoparasites~per~dm^{3}),
             colour = "Species",
             linetype = "Site Elevation",
             shape = 'Species')+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))
F4e

## Fig4f
## Fig4e
F4f <- ggplot(df3, aes(x =(Anelosimus_adults_subadults), y = (klepto_per_capita), colour = factor(species), linetype = factor(site.elevation)))+
        geom_point(aes(shape=factor(species)), size=4) + 
        stat_smooth(method="glm", formula = y ~ x, method.args = list(family=gaussian(link='identity')), se=F, size = 1.5) +
        scale_linetype_manual(values = c("solid", "dashed"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#d3ce3d"))+
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
        labs(x = "# adults and subadults", 
             y = "Kleptoparasites per capita",
             colour = "Species",
             linetype = "Site Elevation",
             shape = 'Species')+
        theme_base()+
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.text=element_text(size=16),
              legend.title = element_text(size=20))+
        theme(legend.key.width = unit(3, 'line'))
F4f

multiplot(F4a, F4b, F4c, F4d, F4e, F4f, cols = 2)
# kleptoDens.per.klepto <- glm(log10(klepto_per_capita) ~ log10(klepto_all) + species.site + log10(klepto_all):species.site, 
#                           data = df3, family = gaussian)
# Anova(kleptoDens.per.klepto, type = "II")
# summary(kleptoDens.per.klepto)
# 
# # Adom400
# F2b.adom400 <- glm(Log10KleptoDens ~ log10(klepto_all), 
#                    data = Adom400)
# Anova(F2b.adom400, type = "III")
# plot(KleptoDens_dm3 ~ log10(klepto_all), data = Adom400)
# 
# # Aexi400
# F2b.aexi400 <- glm(Log10KleptoDens ~ log10(klepto_all), 
#                    data = Aexi400)
# Anova(F2b.aexi400, type = "III")
# plot(Log10KleptoDens ~ log10(klepto_all), data = Aexi400)
# 
# # Aexi1150
# F2b.aexi1150 <- glm(Log10KleptoDens ~ log10(klepto_all), 
#                     data = Aexi1150)
# Anova(F2b.aexi1150, type = "III")
# plot(Log10KleptoDens ~ log10(klepto_all), data = Aexi1150)
# 
# # Aele1150
# F2b.aele1150 <- glm(Log10KleptoDens ~ log10(klepto_all), 
#                     data = Aele1150)
# Anova(F2b.aele1150, type = "III")
# plot(Log10KleptoDens ~ log10(klepto_all), data = Aele1150)
# 
# F4d <- ggplot(df3, aes(x = klepto_all, y = (klepto_per_capita), colour = factor(species.site), linetype = factor(site.elevation)))+
#         geom_point(aes(shape=factor(species))) + 
#         stat_smooth(method="glm", method.args = list(family=gaussian), se=F, size = 1.5) +
#         scale_linetype_manual(values = c("solid", "dashed"))+
#         scale_colour_manual(values = c("#000000", "#666666","#999999", "#CCCCCC"))+
#         scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                       labels = trans_format("log10", math_format(10^.x)))+
#         scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                       labels = trans_format("log10", math_format(10^.x)))+
#         labs(x = "# Kleptoparasites", 
#              y = expression(Kleptoparasites~per~dm^{3}),
#              colour = "species.site",
#              linetype = "Site Elevation")+
#         theme_base()+
#         theme(axis.text = element_text(size = 24),
#               axis.title = element_text(size = 24),
#               legend.text=element_text(size=20),
#               legend.title = element_text(size=24))+
#         theme(legend.key.width = unit(3, 'line'))
# 
# F4d
########### Figure 5 ############
hygieneBeh <- read.csv("hygieneBeh.csv")
cont <- hygieneBeh[,1:15]
cont <- cont[,-c(4:7)]

# remove A and G behaviours, cols 4 & 6
cont <- cont[,-c(4,6)]
View(cont)
cont$Species <- as.character(cont$Species)

rows <- grep("Dom", cont$Nest)
Domingo_400 <- cont[rows,]
which(is.na(Domingo_400$Repair == T)) #22
Domingo_400$Repair[22] <- 0

# want to sum each column 4:9
domSum <- rep(0, times = 10)
for(i in c(4:9)){
        domSum[i] <- sum(Domingo_400[,i])
}
domSum <- domSum[4:9]
domSum <- domSum[1:5]
View(domSum)

rows <- grep("Ex_", cont$Nest)
Eximius_400 <- cont[rows,]
exJsSum <- rep(0, times = 10)
for(i in c(4:9)){
        exJsSum[i] <- sum(Eximius_400[,i])
}
exJsSum <- exJsSum[4:9]
exJsSum <- exJsSum[1:5]

rows <- grep("Exi_", cont$Nest)
Eximius_1150 <- cont[rows,]
exSuSum <- rep(0, times = 10)
for(i in c(4:9)){
        exSuSum[i] <- sum(Eximius_1150[,i])
}
exSuSum <- exSuSum[4:9]
exSuSum <- exSuSum[1:5]


rows <- grep("Ele_VL", cont$Nest)
Elegans_1150 <- cont[rows,]
eleSuSum <- rep(0, times = 10)
for(i in c(4:9)){
        eleSuSum[i] <- sum(Elegans_1150[,i])
}
eleSuSum <- eleSuSum[4:9]
eleSuSum <- eleSuSum[1:5]


rows <- grep("Ele_VB", cont$Nest)
Elegans_1600 <- cont[rows,]
eleCoSum <- rep(0, times = 10)
for(i in c(4:9)){
        eleCoSum[i] <- sum(Elegans_1600[,i])
}
eleCoSum <- eleCoSum[4:9]
eleCoSum <- eleCoSum[1:5]
View(eleCoSum)


rows <- grep("Guac", cont$Nest)
Guac_1600 <- cont[rows,]
guacSum <- rep(0, times = 10)
for(i in c(4:9)){
        guacSum[i] <- sum(Guac_1600[,i])
}
guacSum <- guacSum[4:9]
View(guacSum)
guacSum[1] <- 0
guacSum <- guacSum[1:5]

M <- as.table(rbind(domSum, exJsSum, exSuSum, eleSuSum, eleCoSum, guacSum))
dimnames(M) <- list(Species = c("Domingo 400", "Eximius 400", "Eximius 1150","Elegans 1150", "Elegans 1600", "Guacamayos 1600"),
                    Behavour = c("Selfgroom", "Defecate off nest", "Klepto encounter", 
                                 "Debris removal", "Repair basket"))
M <- matrix(M, 6)

chisq.test(M, simulate.p.value = T)
chisq.test(M)

## Get Fig5 Plots
#proportion of each behaviour
propS <- hygieneBeh$Selfgroom / hygieneBeh$Total
propG <- hygieneBeh$Guard_ES / hygieneBeh$Total
propK <- hygieneBeh$Klepto_chase / hygieneBeh$Total
propD <- hygieneBeh$Debris / hygieneBeh$Total
propP <- hygieneBeh$Defecate / hygieneBeh$Total
propR <- hygieneBeh$Repair / hygieneBeh$Total
Nest <- as.data.frame(hygieneBeh$Nest)
Species <- hygieneBeh$Species
propScans <- cbind(Nest,Species, propS, propK,propD,propR, propP)

# get mean proportion of each behaviour by nest
Selfgroom <- as.numeric(tapply(propScans$propS, propScans$`hygieneBeh$Nest`, mean))
Encounter <- as.numeric(tapply(propScans$propK, propScans$`hygieneBeh$Nest`, mean))
Debris <- as.numeric(tapply(propScans$propD, propScans$`hygieneBeh$Nest`, mean))
Defecate <- as.numeric(tapply(propScans$propP, propScans$`hygieneBeh$Nest`, mean))
Repair <- as.numeric(tapply(propScans$propR, propScans$`hygieneBeh$Nest`, mean))
nest <- as.data.frame(levels(hygieneBeh$Nest))

meanBeh <- cbind(nest, Selfgroom, Encounter, Debris, Defecate, Repair)

testGather <- gather(meanBeh, behaviour, proportion, Selfgroom:Repair)
colnames(testGather)[1] <- "Nest"

## separate out each species by locale ##
# grep() searches for partial strings
## domingo
rows <- grep("Dom", testGather$Nest)
domBeh <- testGather[rows,]
## eximius JS
rows <- grep("Ex_", testGather$Nest)
exBehJS <- testGather[rows,]
## eximius Sumaco
rows <- grep("Exi", testGather$Nest)
exBehSu <- testGather[rows,]
## elegans Sumaco
rows <- grep("Ele_VL", testGather$Nest)
eleBehSu <- testGather[rows,]
## elegans Cocodrilos
rows <- grep("Ele_VB", testGather$Nest)
eleBehCo <- testGather[rows,]
## guacamayos Cocodrilos
rows <- grep("Guac", testGather$Nest)
guacBeh <- testGather[rows,]

### Plot each species separately, then use multiplot function to combine
## Colors: #00b979ff" - eximius, "#554236" - domingo,"#d3ce3d" - elegans, "#f1efa5" - guacamayos
p1 <- ggplot(data = domBeh,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#554236") +
        ggtitle("A. domingo 400") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p1

p2 <- ggplot(data = exBehJS,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#00b979ff") +
        ggtitle("A. eximius 400") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p2

p3 <- ggplot(data = exBehSu,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#00b979ff") +
        ggtitle("A. eximius 1150") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p3

p4 <- ggplot(data = eleBehSu,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#d3ce3d") +
        ggtitle("A. elegans 1150") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p4

p5 <- ggplot(data = eleBehCo,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#d3ce3d") +
        ggtitle("A. elegans 1600") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p5

p6 <- ggplot(data = guacBeh,
             aes(x = behaviour, y = (proportion*100))) + 
        geom_bar(stat = 'summary', fun.y = 'mean', fill = "#f1efa5") +
        ggtitle("A. guacamayos 1600m") +
        expand_limits(y=c(0,1.5))+
        scale_x_discrete(name = "Nest Hygiene Behavior") +
        scale_y_continuous(name = "% Exhibiting Behavior")+
        theme_base()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
p6

multiplot(p5, p4, p1, p6, p3, p2, cols = 2)
############# Figure S1 ###########
klepto.elevation <- glmer(klepto_all ~ nest_elevation + (1|species) + (1|site.elevation), data = df, family = poisson)
Anova(klepto.elevation, type = "III")

FS1 <- ggplot(df, aes(x=nest_elevation, y=klepto_all)) + 
        geom_point() + 
        stat_smooth(method="glm", method.args = list(family="poisson"), se=F, size = 1.5, colour = 'black') +
        labs(x = "Nest elevation (m)", 
             y = "# Kleptoparasites")+
        theme_base()+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24))
FS1


########### Figure S3 ################
# Fig S3a
# policing behaviours #
hygieneBeh <- read.csv("hygieneBeh.csv")
hygieneBeh <- hygieneBeh[,c(1:15)]

# hygieneBeh$allProp <- rowSums(tempdf) / hygieneBeh$Total
hygieneBeh$propK <- hygieneBeh$Klepto_chase / hygieneBeh$Total
hygieneBeh$propD <- hygieneBeh$Debris / hygieneBeh$Total
hygieneBeh$PropP <- hygieneBeh$Defecate / hygieneBeh$Total

# allProp <- (tapply(hygieneBeh$allProp, hygieneBeh$Nest, mean))
K <- (tapply(hygieneBeh$propK, hygieneBeh$Nest, mean))
D <- (tapply(hygieneBeh$propD, hygieneBeh$Nest, mean))
P <- (tapply(hygieneBeh$PropP, hygieneBeh$Nest, mean))

behdf <- as.data.frame(cbind(levels(hygieneBeh$Nest), K, D, P))
# behdf$allProp <- as.numeric.factor(behdf$allProp)
behdf$K <- as.numeric.factor(behdf$K)
behdf$D <- as.numeric.factor(behdf$D)
behdf$P <- as.numeric.factor(behdf$P)

s3a <- merge(df, behdf, by.x = 'nest', by.y = 'V1', all = F)
# K, D, P rows 22:24
tempdf <- s3a[,c(22:24)]
s3a$KDP <- rowSums(tempdf)
s3a$KDP_percent <- s3a$KDP * 100

p <- glm(KDP_percent ~ site.elevation/species, data = s3a, family = gaussian)
Anova(p, type = "II")

ggplot(s3a, aes(x=site.elevation, y=KDP_percent, fill = species))+
        geom_boxplot()+
        theme_base()+
        scale_fill_manual(values = c("#00b979ff", "#554236","#d3ce3d", "#f1efa5"))+
        scale_colour_manual(values = c("#00b979ff", "#554236","#f7bb94ff", "#d3ce3d"))+
        scale_x_discrete(name = "Site elevation") +
        scale_y_continuous(name = "% Performing Behaviors")+
        labs(fill = "Species", colour = "Species")+
        theme(axis.text = element_text(size = 24),
              axis.title = element_text(size = 24),
              legend.text=element_text(size=20),
              legend.title = element_text(size=24))


######### Function: plotting figure array #########
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}

multiplot(F4a, F4b, F4c, cols=1)