#Loading necessary libraries
#If libraries not installed already, run: install.packages("library name")
library(caTools)
library(tidyverse)
library(MASS)
library(VGAM)
library(lme4)
library(pscl)
library(corrr)
library(stargazer)
library(glmulti)
library(sandwich)
library(lmtest)
library(glmmTMB)
library("performance")
library(corrplot)
##############################################################################################################

# Data cleaning
# Make sure that the files are in the current directory

# Read csv file containing the sites for model validation
validation <- read_csv("Validation_1.csv")

# Read csv file containing field data
df <- read_csv("Dataset_257_vphl1.csv")

# Data for model development (exclusive lane)
df1 <- filter(df, `Lane configuration` == "Exclusive")
df1 <- filter(df1, !(Location %in% validation$Location))

# Converting 5 mins counts to vphl
df1[c(8:22, 24:25)] <- lapply(df1[c(8:22, 24:25)], function (x) {x*12})
df1$`Red duration` <- df1$`Red duration`/300

#Adding additional columns for data availability
df1$`Total conflicting thru volume` <- df1$`Conflicting thru volume during red` + 
  df1$`Conflicting thru volume during green`
df1$`Total opposing left-turn volume` <- df1$`Opposing left volume during red` + 
  df1$`Opposing left volume during green`
df1$`Total shadowed left-turn volume` <- df1$`Shadowed left volume during red` + 
  df1$`Shadowed left volume during green`
df1$`Total pedestrian volume (Parallel)` <- df1$`Parallel ped volume during red` + 
  df1$`Parallel ped volume during green`
df1$`Total pedestrian volume (Conflicting)` <-  df1$`Conflicting ped volume during red` + 
  df1$`Conflicting ped volume during green`

# Correlation matrix of quantitative variables
# Outputs an image in the same directory
a <- df1[,c(8,6, 9:23, 87:91)]
a <- as.data.frame(a)
a <- filter(a, is.na(`Red duration`) == FALSE)

tiff("test.tiff", units="in", width=30, height=16, res=300)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(a),  tl.cex = 1.3, type = "upper",
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         cl.cex = 1.1,
         number.cex= 1.3,
         order="original",
         col=col(200), # colour of glyphs
         addCoef.col="black", # colour of coefficients
)

dev.off()

# Descriptive statistics of categorical variables
# Outputs html file with table
a <- df1[,c(8,6, 9:23, 87:91)]
a <- as.data.frame(a)
library(stargazer)
stargazer(a, type = "html", digits=1, out="table3.htm")

# datasets to be used for modeling
c <- df1[c(3, 8, 6, 9:23, 87:91, 36,37,40,41,78,79,81,83:86,4)]
c$`Subject Approach speed limit` <- as.numeric(c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <-  as.numeric(c$`Crossing Approach Speed Limit`)
c$`Subject Approach speed limit` <- ifelse(is.na(c$`Subject Approach speed limit`) == TRUE, 35, c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <- ifelse(is.na(c$`Crossing Approach Speed Limit`) == TRUE, 35, c$`Crossing Approach Speed Limit`)

# Dataset for all models except mixed effect and logistic models
d <- c[c(2:9, 17, 19:34)]
names(d) <- make.names(names(d), unique=TRUE)
dfe <- d

# Dataset for logistic regression
d <- c[c(18, 3:9, 17, 19:34)] # For logistic only
names(d) <- make.names(names(d), unique=TRUE)
dfe_logistic <- d

# Dataset for mixed effect model
d <- c[c(35,2:9, 17, 19:34)] # For mixed only
names(d) <- make.names(names(d), unique=TRUE)
dfe_mixed <- d

###################################### Poisson regression #######################################

# Model form 1
mp_e1 <- glm(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
              Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
              Presence.of.parallel.pedestrian.crosswalk   +
              Receiving.lane + Shadowed.left, family = poisson, data = dfe)
summary(mp_e1)

# Model form 2
mp_e2 <- glm(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
               Total.opposing.left.turn.volume +
               Total.shadowed.left.turn.volume+ Total.right.turn.volume + Total.pedestrian.volume..Conflicting., family = poisson, data = dfe)
summary(mp_e2)

# Model form 3
mp_e3 <- glm(formula = Right.turn.on.red.volume ~ 1 + Total.right.turn.volume, family = poisson, data = dfe)
summary(mp_e3)

# Outputs the model results in a table
stargazer(mp_e1,mp_e2,mp_e3,
          single.row = TRUE, 
          type ="html", 
          report = "vc*", 
          header = FALSE, 
          df=FALSE, 
          digits=2, 
          se = NULL,
          intercept.bottom = FALSE,
          out="models.htm"
)

# Checks for overdispersion
mp_e <- mp_e1
mp_e$deviance/mp_e$df.residual

#Residual Plot
mu.hat <- mp_e$fitted.values
stand.resid <- rstandard(model = mp_e, type = "pearson")  

tiff("test.tiff", units="in", width=13, height=7, res=300)

plot(x = mu.hat, y = stand.resid, xlab = expression(hat(mu)), 
     ylab = "Standardized Pearson residuals", main = "Residual Plot",
     ylim = c(min(c(-3, stand.resid)), max(c(3, stand.resid))), pch=19, cex.lab=1.5, cex.axis=2, cex.main=2, cex.sub=2)
abline(h = c(-3,-2,0,2,3), lty = "dotted", col = "red")

dev.off()

pchisq(mp_e$deviance, df=mp_e$df.residual, lower.tail=FALSE)

# There is evidence of overdispersion

# Excess zeros
tiff("test.tiff", units="in", width=16, height=13, res=300)

ggplot(c, aes(Right.turn.on.red.volume)) + geom_histogram(binwidth = 1) +
  labs(x = "",y="", title="") + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = rel(3)))+    ##Change size of text in title
  theme(axis.title.y = element_text(size = rel(1.4)))+    ##Change size of label text in y-axis label
  theme(axis.title.x = element_text(size = rel(1.4)))+    ##Change size of label text in x-axis label
  theme(axis.text.x = element_text(size = rel(5)))+    ##Change size of tick text in y-axis label
  theme(axis.text.y = element_text(size = rel(5)))
dev.off()

################### Negative binomial ##################################

# Model form 1
mnb_e1 <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                  Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                  Presence.of.parallel.pedestrian.crosswalk   +
                  Receiving.lane + Shadowed.left, data = dfe)

summary(mnb_e1)

# Model form 2
mnb_e2 <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                   Total.opposing.left.turn.volume +
                   Total.shadowed.left.turn.volume+ Total.right.turn.volume + Total.pedestrian.volume..Conflicting., data = dfe)

summary(mnb_e2)

# Model form 3
mnb_e3 <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Total.right.turn.volume, data = dfe)

summary(mnb_e3)

################################## QuasiPoisson #####################################

# Model form 1
mqp_e1 <- glm(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
               Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
               Presence.of.parallel.pedestrian.crosswalk   +
               Receiving.lane + Shadowed.left, family = quasipoisson, data = dfe)

summary(mqp_e1)

# Model form 2
mqp_e2 <- glm(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                Total.opposing.left.turn.volume +
                Total.shadowed.left.turn.volume+ Total.right.turn.volume + Total.pedestrian.volume..Conflicting., family = quasipoisson, data = dfe)

summary(mqp_e2)

# Model form 3
mqp_e3 <- glm(formula = Right.turn.on.red.volume ~ 1 + Total.right.turn.volume, family = quasipoisson, data = dfe)

summary(mqp_e3)

check_overdispersion(mp_e)
check_zeroinflation(mznb_e)
compare_performance(mp_e, mnb_e, mh_e, mznb_e, rank = TRUE)
check_collinearity(mnb_e)

##################################### Hurdle model ######################################

# Model form 1
mh_e1 <- hurdle(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                 Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                 Presence.of.parallel.pedestrian.crosswalk   +
                 Receiving.lane + Shadowed.left| Red.duration + Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mh_e1)

# Model form 2
mh_e2 <- hurdle(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                  Total.opposing.left.turn.volume +
                  Total.shadowed.left.turn.volume+ Total.right.turn.volume + 
                  Total.pedestrian.volume..Conflicting.| Red.duration + Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mh_e2)

# Model form 3
mh_e3 <- hurdle(formula = Right.turn.on.red.volume ~ 1 +
                  Total.right.turn.volume| Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mh_e3)

######################### Zero-inflated negative binomial model #######################

# Model form 1
mznb_e1<- zeroinfl(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                    Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                    Presence.of.parallel.pedestrian.crosswalk   +
                    Receiving.lane + Shadowed.left| Red.duration + Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mznb_e1)

# Model form 2
mznb_e2<- zeroinfl(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                     Total.opposing.left.turn.volume +
                     Total.shadowed.left.turn.volume+ Total.right.turn.volume + 
                     Total.pedestrian.volume..Conflicting.| Red.duration + Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mznb_e2)

# Model form 3
mznb_e3<- zeroinfl(formula = Right.turn.on.red.volume ~ 1 +
                     Total.right.turn.volume| Total.right.turn.volume, data = dfe, dist = "negbin")
summary(mznb_e3)

# Vuong test for excess zeros
vuong(mnb_e, mznb_e) 

################################# Mixed-effect model ########################################

dfe_mixed$State <- str_to_title(dfe_mixed$State)

# Model form 1
mm_e1 <- glmer.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                    Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                    Presence.of.parallel.pedestrian.crosswalk   +
                    Receiving.lane + Shadowed.left +
                    (1|State), data=dfe_mixed)


summary(mm_e1)

# Model form 2
mm_e2 <- glmer.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                    Total.opposing.left.turn.volume +
                    Total.shadowed.left.turn.volume+ Total.right.turn.volume + 
                    Total.pedestrian.volume..Conflicting. +
                   (1|State), data=dfe_mixed)


summary(mm_e2)

# Model form 3
mm_e3 <- glmer.nb(formula = Right.turn.on.red.volume ~ 1 + Total.right.turn.volume +
                    (1|State), data=dfe_mixed)


summary(mm_e3)


# Random effect plot
plot_model(mm_e3, type = "re") +theme_sjplot2(base_size = 30, base_family = "")

##################################### Logistic regression ########################################

dfe_logistic$Right.turn.on.red.percent <- dfe_logistic$Right.turn.on.red.percent/100

# Model form 1
ml_e1 <- glm(Right.turn.on.red.percent ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
               Shadowed.left.volume.during.red + Conflicting.ped.volume.during.red +
               Presence.of.parallel.pedestrian.crosswalk   +
               Receiving.lane + Shadowed.left,family = binomial, weights = Total.right.turn.volume,
            data = dfe_logistic)

summary(ml_e1)

# Model form 2
ml_e2 <- glm(Right.turn.on.red.percent ~ 1 + Red.duration + Total.conflicting.thru.volume + 
               Total.opposing.left.turn.volume +
               Total.shadowed.left.turn.volume+ 
               Total.pedestrian.volume..Conflicting.,family = binomial, weights = Total.right.turn.volume,
             data = dfe_logistic)

summary(ml_e2)

# Model form 3
ml_e3 <- glm(Right.turn.on.red.percent ~ 1 + Red.duration, family = binomial, weights = Total.right.turn.volume,
             data = dfe_logistic)

summary(ml_e3)

######################### Separate model for Interchange/not? #################################################

# Analysis using data for interchanges only
# Same steps for data preparation as before
validation <- read_csv("Validation_1.csv")
df <- read_csv("Dataset_257_vphl1.csv")

# Data for interchanges with RTOR from exclusive RT lane
df1 <- filter(df, Type == "Interchange", `Lane configuration` == "Exclusive")
df1 <- filter(df1, !(Location %in% validation$Location))

df1[c(8:22, 24:25)] <- lapply(df1[c(8:22, 24:25)], function (x) {x*12})
df1$`Red duration` <- df1$`Red duration`/300
df1$`Total conflicting thru volume` <- df1$`Conflicting thru volume during red` + 
  df1$`Conflicting thru volume during green`
df1$`Total opposing left-turn volume` <- df1$`Opposing left volume during red` + 
  df1$`Opposing left volume during green`
df1$`Total shadowed left-turn volume` <- df1$`Shadowed left volume during red` + 
  df1$`Shadowed left volume during green`
df1$`Total pedestrian volume (Parallel)` <- df1$`Parallel ped volume during red` + 
  df1$`Parallel ped volume during green`
df1$`Total pedestrian volume (Conflicting)` <-  df1$`Conflicting ped volume during red` + 
  df1$`Conflicting ped volume during green`

c <- df1[c(3, 8, 6, 9:23, 87:91, 36,37,40,41,78,79,81,83:86,4)]
#c <- filter(c, is.na(`Red duration`) == FALSE)
c$`Subject Approach speed limit` <- as.numeric(c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <-  as.numeric(c$`Crossing Approach Speed Limit`)
c$`Subject Approach speed limit` <- ifelse(is.na(c$`Subject Approach speed limit`) == TRUE, 35, c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <- ifelse(is.na(c$`Crossing Approach Speed Limit`) == TRUE, 35, c$`Crossing Approach Speed Limit`)

names(c)
d <- c[c(2:9, 17, 19:34)]
names(d) <- make.names(names(d), unique=TRUE)
dfe_interchange <- d
d <- c[c(18, 3:9, 17, 19:34)] # For logistic only
names(d) <- make.names(names(d), unique=TRUE)
dfe_logistic <- d

# Model form 1
mznb_interchange1 <- zeroinfl(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                                Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                                Presence.of.parallel.pedestrian.crosswalk   +
                                Receiving.lane + Shadowed.left| Red.duration + Total.right.turn.volume, data = dfe_interchange, dist = "negbin")

summary(mznb_interchange1)

# Model form 2
mnb_interchange2 <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                             Total.opposing.left.turn.volume +
                             Total.shadowed.left.turn.volume+ Total.right.turn.volume +
                             Total.pedestrian.volume..Conflicting., data = dfe_interchange)


summary(mnb_interchange2)

# Model form 3
dfe_logistic$Right.turn.on.red.percent <- dfe_logistic$Right.turn.on.red.percent/100
ml_interchange3 <- glm(Right.turn.on.red.percent ~ 1 + Red.duration, family = binomial, weights = Total.right.turn.volume,
                        data = dfe_logistic)

summary(ml_interchange3)


# Analysis using data for both intersections and interchanges
validation <- read_csv("Validation_1.csv")
df <- read_csv("Dataset_257_vphl1.csv")
df1 <- filter(df, `Lane configuration` == "Exclusive")
df1 <- filter(df1, !(Location %in% validation$Location))
sites <- df1 %>% group_by(Location) %>% summarise(n= n())

df1[c(8:22, 24:25)] <- lapply(df1[c(8:22, 24:25)], function (x) {x*12})
df1$`Red duration` <- df1$`Red duration`/300
df1$`Total conflicting thru volume` <- df1$`Conflicting thru volume during red` + 
  df1$`Conflicting thru volume during green`
df1$`Total opposing left-turn volume` <- df1$`Opposing left volume during red` + 
  df1$`Opposing left volume during green`
df1$`Total shadowed left-turn volume` <- df1$`Shadowed left volume during red` + 
  df1$`Shadowed left volume during green`
df1$`Total pedestrian volume (Parallel)` <- df1$`Parallel ped volume during red` + 
  df1$`Parallel ped volume during green`
df1$`Total pedestrian volume (Conflicting)` <-  df1$`Conflicting ped volume during red` + 
  df1$`Conflicting ped volume during green`

names(df1)
c <- df1[c(3, 8, 6, 9:23, 87:91, 36,37,40,41,78,79,81,83:86,4, 80)]
#c <- filter(c, is.na(`Red duration`) == FALSE)
c$`Subject Approach speed limit` <- as.numeric(c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <-  as.numeric(c$`Crossing Approach Speed Limit`)
c$`Subject Approach speed limit` <- ifelse(is.na(c$`Subject Approach speed limit`) == TRUE, 35, c$`Subject Approach speed limit`)
c$`Crossing Approach Speed Limit` <- ifelse(is.na(c$`Crossing Approach Speed Limit`) == TRUE, 35, c$`Crossing Approach Speed Limit`)

names(c)
d <- c[c(2:9, 17, 19:34, 36)]
names(d) <- make.names(names(d), unique=TRUE)
d$Type <- as.factor(d$Type)
d$Type <- relevel(d$Type, ref = "Intersection")
dfe_both <- d
d <- c[c(18, 3:9, 17, 19:34,36)] # For logistic only
names(d) <- make.names(names(d), unique=TRUE)
d$Type <- as.factor(d$Type)
d$Type <- relevel(d$Type, ref = "Intersection")
dfe_logistic <- d

# Model form 1
# without categorical variable
mznb <- zeroinfl(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Conflicting.thru.volume.during.red  +
                                Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                                Presence.of.parallel.pedestrian.crosswalk   +
                                Receiving.lane + Shadowed.left| Red.duration + Total.right.turn.volume, data = dfe_both, dist = "negbin")
# with categorical variable
mznb_type <- zeroinfl(formula = Right.turn.on.red.volume ~ 1 + Type + Red.duration + Conflicting.thru.volume.during.red  +
                   Shadowed.left.volume.during.red+ Total.right.turn.volume + Conflicting.ped.volume.during.red +
                   Presence.of.parallel.pedestrian.crosswalk   +
                   Receiving.lane + Shadowed.left| Red.duration + Total.right.turn.volume, data = dfe_both, dist = "negbin")

summary(mznb)

# Model form 2
# without categorical variable
mnb <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Red.duration + Total.conflicting.thru.volume + 
                             Total.opposing.left.turn.volume +
                             Total.shadowed.left.turn.volume+ Total.right.turn.volume +
                Total.pedestrian.volume..Conflicting., data = dfe_both)

# with categorical variable
mnb_type <- glm.nb(formula = Right.turn.on.red.volume ~ 1 + Type+ Red.duration + Total.conflicting.thru.volume + 
                Total.opposing.left.turn.volume +
                Total.shadowed.left.turn.volume+ Total.right.turn.volume +
                Total.pedestrian.volume..Conflicting., data = dfe_both)

summary(mnb)

# Model form 3
dfe_logistic$Right.turn.on.red.percent <- dfe_logistic$Right.turn.on.red.percent/100

# without categorical variable
ml <- glm(Right.turn.on.red.percent ~ 1 + Red.duration, family = binomial, weights = Total.right.turn.volume,
                       data = dfe_logistic)

# with categorical variable
ml_type <- glm(Right.turn.on.red.percent ~ 1 + Type + Red.duration, family = binomial, weights = Total.right.turn.volume,
          data = dfe_logistic)

summary(ml)
