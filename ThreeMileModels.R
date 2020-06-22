# Three Mile Manuscript

# Get and load packages
loadrequired <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, library, character.only = TRUE)
}

packages <- c("readxl", "ggplot2", "dplyr", "nicheROVER", "viridis", "questionr", "MixSIAR", "ggpubr", 
              "cowplot", "tictoc")
loadrequired(packages)

isodata <- read_xlsx("All fish 2012 Isotope analysis.xlsx")

# Clean up the grouping designations

isodata$Class[isodata$Class == "< 40"] <- "<40"
isodata$Class[isodata$Class == "s"] <- "S"
isodata$Time[isodata$Time == "FALL"] <- "Fall"

# Define predators and prey
isodata$Type <- ifelse(isodata$Class == "M" | isodata$Class == "P" | isodata$Class == "Q" | isodata$Class == "S" | 
                          isodata$Class == "T", "Predator", "Prey")

# Instrumental + analytical error reported by SIPERG lab C±0.07 and N±0.11

# Quick exploratory plot
ggplot(data = isodata) + geom_errorbar(aes(x=Carbon, ymax=Nitrogen+0.11, ymin=Nitrogen-0.11), color="gray30") +
                  geom_errorbarh(aes(y=Nitrogen, xmin=Carbon-0.07, xmax=Carbon+0.07), color="gray30")+
  geom_point(aes(x=Carbon, y=Nitrogen, shape=Time, color=Class), size=3, alpha=0.7) + 
  facet_grid(~Type+Time) + scale_shape_manual(values=c(19,17)) +
  theme_bw(base_size = 12) + theme(strip.background = element_rect(colour="black", fill="white")) +
  ylab(expression(paste(delta^15, "N (\u2030)",sep=""))) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep=""))) +
  scale_color_viridis_d(option = "D")

#ggsave("biplot.png", dpi=400, scale=1.25)


# Get just predators
pred.spring <- isodata[isodata$Type %in% "Predator" & isodata$Time %in% "Spring",]
pred.fall <- isodata[isodata$Type %in% "Predator" & isodata$Time %in% "Fall",]

###########################################################################################################
# Models
# stock-quality (S-Q), quality-preferred (Q-P), preferred-memorable (P-M), and memorable-trophy (M-T)
isodata$Class <- factor(isodata$Class, levels = c("<40", "40-59", "60-79", "80-99", "100-120", "S", "Q", "P", "M", "T"))
rawlength <- read.csv("threemilelength.csv")

summary(glm(Nitrogen ~ Species*mm, data=rawlength))
summary(glm(Nitrogen ~ Species*Class, data=isodata))

ggplot(data = isodata) + geom_point(aes(x=Class, y=Nitrogen, color=Species)) + geom_smooth(aes(x=Class, y=Nitrogen, group=Species, color=Species))

d15npanel <- ggplot(data = isodata) + geom_boxplot(aes(x=factor(Class, 
             level=c("<40", "40-59", "60-79", "80-99", "100-120", "S", "Q", "P", "M", "T")), y=Nitrogen, fill=Time, color=Time)) +   facet_grid(~Species) + 
  theme_classic() +
 theme(panel.border = element_rect(color="black", fill=NA))+
 theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  ylab(expression(paste(delta^15, "N (\u2030)",sep=""))) + xlab("Length class") +
  scale_fill_manual(values = c("lightblue", "gray90")) +
  scale_color_manual(values = c("grey60", "black")) + 
theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  plot.margin = margin(1, 2, 2, 2, "pt"))


d13Cpanel <- ggplot(data = isodata) + geom_boxplot(aes(x=factor(Class, 
             level=c("<40", "40-59", "60-79", "80-99", "100-120", "S", "Q", "P", "M", "T")), y=Carbon, fill=Time, color=Time)) + 
  facet_grid(~Species) + 
  theme_classic() +
 theme(panel.border = element_rect(color="black", fill=NA))+
 theme(axis.text.x = element_blank(),
  plot.margin = margin(1, 2, 0, 2, "pt")) +
  ylab(expression(paste(delta^13, "C (\u2030)",sep=""))) + xlab(NULL) +
  scale_fill_manual(values = c("lightblue", "gray90")) +
  scale_color_manual(values = c("grey60", "black"))
  
isoplot <- plot_grid(d13Cpanel + theme(legend.position = "none"),
                     d15npanel + theme(legend.position = "none"), 
          align = "v", 
          axis = "b", 
          nrow=2,
          rel_heights = c(7/16, 9/16),
                          labels = c("A", "B"))
legend <- get_legend(
  # create some space to the left of the legend
  d15npanel + theme(legend.box.margin = margin(0, 0, 0, 2))
)
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
isopanelplot <- plot_grid(isoplot, legend, ncol = 2, rel_widths = c(7/8, 1/8))
save_plot("SI_Figure_1.pdf", isopanelplot, dpi=1000, base_height = 7.5)
#ggsave("d15Npanel.png", dpi=400, width = 8.5)

ggplot(isodata) + geom_boxplot(aes(x=Class, y=Carbon)) + 
  facet_grid(Time~Species) + 
  theme_classic() +
 theme(panel.border = element_rect(color="black", fill=NA))+
 theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  ylab(expression(paste(delta^13, "C (\u2030)",sep=""))) + xlab("Size Class")
  

# ggsave("d13Cpanel.png", dpi=400, width = 8.5)

###########################################################################################################

# Setup for Swanson models:SPRING
nsamples <- 1000
system.time({
    iso.par <- tapply(1:nrow(pred.spring), pred.spring$Species, function(ii) niw.post(nsamples = nsamples, 
        X = pred.spring[ii, 2:3]))
})

summary.factor(pred.spring$Species) # How many colors do we need?
clrs <- viridis(6, option = "inferno", end=0.8)  # colors for each species
clrs2 <- c("gray25", "gray30", "gray35", "gray40", "gray45", "gray50")
# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, 0.5, 0.1) + 0.1, mfrow = c(1, 3))
niche.par.plot(iso.par, col = clrs, plot.index = 1)
niche.par.plot(iso.par, col = clrs, plot.index = 2)
niche.par.plot(iso.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(iso.par), fill = clrs)

# all mu (del15N, del13C)
niche.par.plot(iso.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topright", legend = names(iso.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1) + 0.1)
niche.par.plot(iso.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(iso.par), fill = clrs)

# 2-d projections of 10 niche regions
nsamples <- 10
iso.par <- tapply(1:nrow(pred.spring), pred.spring$Species, function(ii) niw.post(nsamples = nsamples, 
    X = pred.spring[ii, 2:3]))

# format data for plotting
iso.data <- tapply(1:nrow(pred.spring), pred.spring$Species, function(ii) X = pred.spring[ii, 2:3])
iso.data <- lapply(iso.data, as.data.frame)
niche.plot(niche.par = iso.par, niche.data = iso.data, pfrac = 0.05, 
           iso.names = expression(delta^{13} * C, delta^{15} * N), col = clrs, xlab = expression(paste("Isotope Ratio (\u2030)")))

nsamples <- 1000
# Overlap calculation
over.stat <- overlap(iso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
    0.99))

# The mean overlap metrics calculated across iteratations for both niche
# region sizes (alpha = .95 and alpha = .99) 
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
round(over.mean, 2)

over.cred <- apply(over.stat * 100, c(1:2, 4), quantile, prob = c(0.025, 0.975), 
    na.rm = TRUE)
round(over.cred[, , , 1])  # display alpha = .95 niche region

# Overlap plot at an alpha of 0.95
over.stat <- overlap(iso.par, nreps = nsamples, nprob = 1000, alpha = 0.95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "cyan2", equal.axis = TRUE, 
    xlab = "Spring Overlap Probability (%): Niche Region Size: 95%")

# Only HSB collected in predator size in fall

#####################################################################################################################
# 3.4 N
# C 1
# MixSIAR models
inverts <- read_excel("Zoop.xlsx")
inverts$Type <- "Prey"

allprey <- rbind(isodata[isodata$Type == "Prey", ], inverts)

allprey$Source <- paste(allprey$Species, allprey$Class, allprey$Time)
allprey$Source <- recode(allprey$Source, `ZOO Zoo` = "ZOO")
ggplot(allprey) + geom_point(aes(x=Carbon, y=Nitrogen, color=Species, shape=Time)) + 
  scale_color_viridis_d()
allprey$Source <- paste(allprey$Species, allprey$Class, allprey$Time)
allprey <- allprey[,c(8,1:7)]
write.csv(allprey, file="allprey.csv", row.names = FALSE)
hsb.spring <- isodata[isodata$Species == "HSB" & isodata$Type == "Predator" & isodata$Time == "Spring", ]
hsb.predators <- isodata[isodata$Species == "HSB" & isodata$Type == "Predator", ]
write.csv(hsb.predators, "hsbonlypredators.csv", row.names = FALSE)

prey.spring <- allprey[allprey$Time == "Spring", ]
prey.spring <- prey.spring[!(prey.spring$Identifier=="188"),] # Drop lone observation

#prey.spring <- as.data.frame(prey.spring)
prey.spring$Source <- paste(prey.spring$Species,prey.spring$Class)

prey.spring <- prey.spring[ ,c(8, 2:3,5)]
prey.spring$Source <- recode(prey.spring$Source, `ZOO Zoo` = "ZOO")
prey.spring$Species <- as.factor(prey.spring$Species)
prey.spring$Class <- as.factor(prey.spring$Class)
prey.spring$Source <- as.factor(prey.spring$Source)

prey.spring <- prey.spring %>% 
  group_by(Species) %>% 
  summarise(MeanCarbon = mean(Carbon), MeanNitrogen= mean(Nitrogen), SDGCarbon = sd(Carbon), SDNitrogen = sd(Nitrogen), n=length(Species))

write.csv(prey.spring, file="allpreyspring.csv", row.names = FALSE)

discrimspring <- unique(prey.spring[,8])
discrimall <- unique(allprey[,c(4:6)])
discrimall$Source <- paste(discrimall$Species, discrimall$Class, discrimall$Time)
discrimall <- discrimall[,c(4)]
discrimall$MeanCarbon <- 0.8
discrimall$SDCarbon <- 1.1
discrimall$MeanNitrogen <- 3.4
discrimall$SDNitrogen <- 0.5
write.csv(discrimall, file="discrimall.csv", row.names = FALSE)
#names(discrimspring) <- "Species"
#discrimspring$Source <- recode(discrimspring$Source, `ZOO Zoo` = "ZOO")
discrimspring$MeanCarbon <- 0.8
discrimspring$SDCarbon <- 1.1
discrimspring$MeanNitrogen <- 3.4
discrimspring$SDNitrogen <- 0.5

write.csv(discrimspring, file="discrimspring.csv", row.names = FALSE)

mixhsbspring <- load_mix_data(filename="hsbonly.csv",
					 iso_names=c("Carbon","Nitrogen"),
					 factors=NULL,
					 fac_random=NULL,
					 fac_nested=NULL,
					 cont_effects=NULL)

# Load source data
source <- load_source_data(filename="allprey.csv",
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="raw",
						   mixhsbspring)

discr <- load_discr_data(filename="discrimall.csv", mixhsbspring)

calc_area(source=source,mix=mixhsbspring, discr = discr)
plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mixhsbspring, source)

tic()
hsbspring.jags <- run_model(run="normal", mixhsbspring, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)
toc()
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)


plot_data_two_iso(hsbspring.jags, mixhsbspring, source)
output_JAGS(hsbspring.jags, mixhsbspring, source, output_options)

####################################################################################################

# Tests for independence
chisq.test(xtabs(~ Time + Class, data=isodata), simulate.p.value=TRUE)
chisq.test(xtabs(~ Carbon + Nitrogen, data=isodata), simulate.p.value=TRUE)
chisq.test(xtabs(~ Species + Class, data=isodata), simulate.p.value=TRUE)

# Count n by species and class
summarytable <- isodata %>%
 count(Species, Class, sort = F)

# Tests for normality
ggdensity(isodata$Carbon, 
          main = "Density plot of d13C",
          xlab = "d13C")

ggdensity(isodata$Nitrogen, 
          main = "Density plot of d15N",
          xlab = "d15N")

###############################################################################################################

# Mixing models for hsb
hsb.mix <- as.matrix(isodata[isodata$Species == "HSB" & isodata$Type == "Predator", 2:3]) 

colnames(hsb.mix) <- c("d13C", "d15N")

source.means <- isodata[isodata$Type == "Prey", c(2:5)] %>%
  group_by(Species, Class) %>% # This line then groups the data by ID
  summarize_all(.funs = c(mean="mean")) # This line then calculates the mean d15C and d15N by the groups we created above.

source.means2 <- as.matrix(source.means[-c(1,4), 3:4])

source.names <- source.means[-c(1,4), ] %>%
  unite(Prey, Species, Class, remove = FALSE)

source.names2 <- as.vector(source.names$Prey)

# SDS
source.sds <- isodata[isodata$Type == "Prey", c(2:5)] %>%
  group_by(Species, Class) %>% # This line then groups the data by ID
  summarize_all(.funs = c(sd="sd"))

source.sds2 <- as.matrix(source.sds[-c(1,4), 3:4])

hsb.simmr.in = simmr_load(mixtures=hsb.mix,
                     source_names=source.names2,
                     source_means=source.means2,
                     source_sds=source.sds2)

plot(hsb.simmr.in, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), # \u2030 is the unicode for the permil symbol
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), # This gives it a formatted y-axis label
     title="Isospace plot of potential hsb food sources") # This gives it a plot title

# Note: The point or points labeled "Mixtures" are the isotopic signatures of the organism for which you are modeling diet. The Mixture in our case is the trout.

# The command simmr_mcmc creates an isotope mixing model using Markov chain Monte Carlo.
hsb.simmr.out = simmr_mcmc(hsb.simmr.in)

summary(hsb.simmr.out,type='diagnostics')


# This command creates a box and whisker plot of the potential food sources
compare_sources(hsb.simmr.out)

###################################################################################
# Mixing models for yeb

yeb.mix <- as.matrix(spring[spring$Species == "YEB", 2:3]) 

colnames(yeb.mix) <- c("d13C", "d15N")

yeb.simmr.in = simmr_load(mixtures=yeb.mix,
                     source_names=source.names,
                     source_means=source.means2,
                     source_sds=source.sds2)



plot(yeb.simmr.in, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), # \u2030 is the unicode for the permil symbol
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), # This gives it a formatted y-axis label
     title="Isospace plot of potential YEB food sources") # This gives it a plot title


# The command simmr_mcmc creates an isotope mixing model using Markov chain Monte Carlo.
yeb.simmr.out = simmr_mcmc(yeb.simmr.in)

summary(yeb.simmr.out,type='diagnostics')


# This command creates a box and whisker plot of the potential food sources
compare_sources(yeb.simmr.out)

#####################################################
# Test for differences between size classes

kruskal.test(Nitrogen ~ Source, data=allprey)
kruskal.test(Carbon ~ Source, data=allprey)
library(FSA)
# Order groups by median


dunnNitrogen <- dunnTest(Nitrogen ~ Source,
              data=allprey,
              method="bh") 

dunnCarbon <- dunnTest(Carbon ~ Source,
              data=allprey,
              method="bh")
write.csv(file="dunnnitrogen.csv", dunnNitrogen[2])
write.csv(file="dunncarbon.csv", dunnCarbon[2])

# Combining sources for which there are no significant differences between size classes across time within the same species
allprey$Source <- ifelse(allprey$Species == "BLC", "BLC <40-120", 
                 ifelse(allprey$Species == "BLG" & allprey$Class != "<40", "BLG 40-120", 
                        ifelse(allprey$Species == "BLG" & allprey$Class == "<40", "BLG <40",
                        ifelse(allprey$Species == "LMB", "LMB 40-120",
                               ifelse(allprey$Species == "WHC", "WHC 60-120",
                                      ifelse(allprey$Species == "ZOO", "ZOO", 
                                             ifelse(allprey$Species == "YEB" & allprey$Class != "100-120", "YEB 40-79", 
                                                    ifelse(allprey$Species == "YEB" & allprey$Class == "100-120", "YEB 100-120"                                                              ,ifelse(allprey$Species == "BEN", "BEN", paste())))))))))

# New dual isotope plot
#allprey <- allprey[!(allprey$Identifier=="188"),] # Drop lone observation
write.csv(allprey, "allprey.csv", row.names=FALSE)
pal <- c("#FDE725FF", "#46337EFF", "#440154FF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FCA50AFF")

# Figure 1
ggplot(data = allprey) + 
  geom_errorbar(aes(x=Carbon, ymax=Nitrogen+0.11, ymin=Nitrogen-0.11), color="grey50") +
  geom_errorbarh(aes(y=Nitrogen, xmin=Carbon-0.07, xmax=Carbon+0.07), color="grey50") +
  geom_point(aes(x=Carbon, y=Nitrogen, color=Source, shape=Source), size=3) + 
  # Order: BEN, BLC, BLG <40, BLG 40-120, LMB, WHC, YEB 100-120, YEB 40-79, ZOO
  scale_shape_manual(values=c(16,0,17,17,18,15,8,8,9), name="Class",
                       labels=c("BLC < 40-120 mm","BLG < 40 mm", "BLG 40-120 mm","LMB 40-120 mm", "WHC 60-120 mm", "YEB 40-79 mm",
                                "YEB 100-120 mm", "BEN", "ZOO"),
                       breaks=c("BLC <40-120", "BLG <40", "BLG 40-120","LMB 40-120", "WHC 60-120", "YEB 40-79", "YEB 100-120", "BEN", "ZOO")) +
  theme_bw(base_size = 12) +
  ylab(expression(paste(delta^15, "N (\u2030)",sep=""))) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep=""))) +
  scale_color_manual(values=pal,
   name="Class",
                       labels=c("BLC < 40-120 mm","BLG < 40 mm", "BLG 40-120 mm","LMB 40-120 mm", "WHC 60-120 mm", "YEB 40-79 mm",
                                "YEB 100-120 mm", "BEN", "ZOO"),
                       breaks=c("BLC <40-120", "BLG <40", "BLG 40-120","LMB 40-120", "WHC 60-120", "YEB 40-79", "YEB 100-120", "BEN", "ZOO"))
ggsave("Figure_1.pdf", dpi=1000, height = 6, width = 6)

# Observations by class after pooling
summary.factor(allprey$Source)

allprey %>%
  group_by(Source) %>%
  summarize(mean_d13C = mean(Carbon), mean_d15N = mean(Nitrogen), d13C_min = min(Carbon), d15N_min = min(Nitrogen),
            d13C_max = max(Carbon), d15N_max = max(Nitrogen))


discrimall <- unique(allprey[,1])
discrimall$MeanCarbon <- 0.8
discrimall$SDCarbon <- 1.1
discrimall$MeanNitrogen <- 3.4
discrimall$SDNitrogen <- 0.5
write.csv(discrimall, file="discrimall.csv", row.names = FALSE)

######################################################################
# Mixing models for hsb

allprey <- read.csv("allprey.csv")

hsbT.mix <- as.matrix(isodata[isodata$Species == "HSB" & isodata$Class == "T", 2:3]) 

colnames(hsbT.mix) <- c("d13C", "d15N")

source.means <- allprey[, c(8,2,3)] %>%
  group_by(Source) %>% # This line then groups the data by ID
  summarize_all(.funs = c(mean="mean")) # This line then calculates the mean d15C and d15N by the groups we created above.

source.means2 <- as.matrix(source.means[, 2:3])

source.names <- source.means[ , 1]

source.names2 <- as.vector(source.names$Source)

# SDS
source.sds <- allprey[, c(8,2,3)] %>%
  group_by(Source) %>% # This line then groups the data by ID
  summarize_all(.funs = c(sd="sd"))

source.sds2 <- as.matrix(source.sds[, 2:3])

hsbT.simmr.in = simmr_load(mixtures=hsbT.mix,
                     source_names=source.names2,
                     source_means=source.means2,
                     source_sds=source.sds2)

plot(hsbT.simmr.in, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), # \u2030 is the unicode for the permil symbol
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), # This gives it a formatted y-axis label
     title="Isospace plot of potential hsb food sources") # This gives it a plot title

# Note: The point or points labeled "Mixtures" are the isotopic signatures of the organism for which you are modeling diet. The Mixture in our case is the trout.

# The command simmr_mcmc creates an isotope mixing model using Markov chain Monte Carlo.
hsbT.simmr.out = simmr_mcmc(hsbT.simmr.in)

summary(hsbT.simmr.out,type='diagnostics')


# This command creates a box and whisker plot of the potential food sources
compare_sources(hsbT.simmr.out)

