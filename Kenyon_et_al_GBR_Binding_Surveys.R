# Kenyon et al. (2024) Binding Surveys GBR
# Accepted by GCB

# Load libraries 
library(tidyverse)
library(emmeans)
library(glmmTMB)
library(gridExtra)
library(MASS)
library(DHARMa)
library(cowplot)
library(MuMIn) 
library(performance) # for detecting collinearity between predictor variables
library(see)
library(ordinal) # cumulative link models
library(flextable) # for saving tables
library(officer) # for saving tables

# Writing a function for saving the tables:

save_tables_to_docx <- function(emmeans_obj, 
                                emmeans_filename, 
                                contrasts_filename) {
  # Extract emmeans and contrasts
  emmeans_df <- emmeans_obj$emmeans %>% as.data.frame()
  contrasts_df <- emmeans_obj$contrasts %>% as.data.frame()
  
  # Set flextable defaults
  set_flextable_defaults(font.family = "Times New Roman", font.size = 9,
                         padding.bottom = 0, padding.top = 0)
  
  # Create flextables
  emmeans_table <- flextable(emmeans_df) %>% 
    colformat_double(big.mark = ",", digits = 3, na_str = "N/A")
  
  contrasts_table <- flextable(contrasts_df) %>% 
    colformat_double(big.mark = ",", digits = 3, na_str = "N/A")
  
  # Save flextables to Word documents
  save_as_docx(emmeans_table, path = emmeans_filename)
  save_as_docx(contrasts_table, path = contrasts_filename)
}


# Load data 

dat <- read.csv("Binding_Survey_Data_v8.csv", header =TRUE)

# Data formatting ----

dat2 <- dat %>% dplyr::select("Trip", "Location", "Exposure", "Depth", "Site", "Date", 
                              "Quadrat2", "Number_coral_recruits_less5cm","Survey_depth_m", 
                              "Quadrat_angle2", "Survey_layer",
                              "Patch_width_m", "Patch_length_m", "Patch_area_m2", 
                              "distance_centre_to_closest_patch_edge_m_2", 
                              "Average_bed_depth_cm", "bed_depth_1", "bed_depth_2",                       
                              "bed_depth_3","bed_depth_4", "bed_depth_5", "bed_depth_6", 
                              "bed_depth_7", "bed_depth_8", "Piece",
                              "Unique_piece", "Layer_in_rubble_bed","Percent_hidden_in_bed", 
                              "Number_contact_points2","Max_gap_size_to_right_cm2", "Stable_0.1_v2",
                              "Stable_0.2", "Stable_description2",
                              "Widest_span_cm2", "Morphology2", "Num_branches2","Num_branches3_0_for_other",
                              "Interlocking_morphology_or_not.Branched_Plate", 
                              "Number_binds2", "Bound", "Bound_to_how_many",
                              "Binder1.broad_cat", "Binder1_length_cm",
                              "Binder2.broad_cat", "Binder2_length_cm", 
                              "Binder3.broad_cat", "Binder3_length_cm",
                              "Binder4.broad_cat", "Binder4_length_cm",
                              "Binder5.broad_cat", "Binder5_length_cm",
                              "Binder6.broad_cat", "Binder6_length_cm",
                              "Binder7.broad_cat", "Binder7_length_cm",
                              "Binder8.broad_cat", "Binder8_length_cm",
                              "Binder9.broad_cat", "Binder9_length_cm",
                              "Binder10.broad_cat", "Binder10_length_cm",
                              "Binder11.broad_cat", "Binder11_length_cm",
                              "Binder12.broad_cat", "Binder12_length_cm",
                              "Binder13.broad_cat", "Binder13_length_cm",
                              "Binder14.broad_cat", "Binder14_length_cm",
                              "Binder15.broad_cat", "Binder15_length_cm")

cols1 <- c("Trip", "Location", "Exposure", "Depth", "Site", 
           "Quadrat2","Quadrat_angle2", "Survey_layer",
           "Piece", "Unique_piece", "Layer_in_rubble_bed",
           "Percent_hidden_in_bed", "Stable_description2",
           "Morphology2", "Unique_piece",
           "Interlocking_morphology_or_not.Branched_Plate", 
           "Binder1.broad_cat", 
           "Binder2.broad_cat", 
           "Binder3.broad_cat", 
           "Binder4.broad_cat", 
           "Binder5.broad_cat", 
           "Binder6.broad_cat", 
           "Binder7.broad_cat", 
           "Binder8.broad_cat", 
           "Binder9.broad_cat",
           "Binder10.broad_cat", 
           "Binder11.broad_cat", 
           "Binder12.broad_cat", 
           "Binder13.broad_cat", 
           "Binder14.broad_cat", 
           "Binder15.broad_cat") 

dat2 <- dat2 %>% mutate_at(cols1, factor)

# Create an exposure depth variable

dat2$exposuredepth <- paste(dat2$Exposure, dat2$Depth, sep = "_")
dat2$exposuredepth <- as.factor(dat2$exposuredepth)
levels(dat2$exposuredepth)

# Create a location exposure depth variable

dat2$locexpdepth <- paste(dat2$Location, dat2$exposuredepth, sep = "_")
dat2$locexpdepth <- as.factor(dat2$locexpdepth)
levels(dat2$locexpdepth)

# Format levels

dat2$locexpdepth <- factor(dat2$locexpdepth, levels = c("Keppels (inshore)_Sheltered_Shallow",
                                                        "Keppels (inshore)_Exposed_Shallow",
                                                        "Heron (offshore)_Sheltered_Reef flat",
                                                        "Heron (offshore)_Sheltered_Shallow",
                                                        "Heron (offshore)_Sheltered_Deep",
                                                        "Heron (offshore)_Intermediate_Shallow",
                                                        "Heron (offshore)_Intermediate_Deep",
                                                        "Heron (offshore)_Exposed_Shallow",
                                                        "Heron (offshore)_Exposed_Deep",
                                                        "Heron (offshore)_Intermediate_BP"))

dat2$Site <- factor(dat2$Site, levels = c("Halfway Sheltered", "Humpy Sheltered", "Clam Bay",
                                          "Halfway Exposed", "Humpy Exposed", "Red Beach",
                                          "Coral Gardens", "Halfway", "First Point",
                                          "Eco_1", "Eco_2", "Blue Pools"))

dat2$Location <- factor(dat2$Location, levels = c("Keppels (inshore)",
                                                  "Heron (offshore)"),
                        labels = c("Keppels (inshore)", 
                                   "Heron (offshore)"))

dat2$sitedepth <- paste(dat2$Site, dat2$Depth, sep = "_")
dat2$sitedepth <- as.factor(dat2$sitedepth)
levels(dat2$sitedepth)

dat2$sitedepth <- factor(dat2$sitedepth, levels = c("Humpy Sheltered_Shallow",
                                                    "Halfway Sheltered_Shallow",
                                                    "Clam Bay_Shallow",
                                                    "Humpy Exposed_Shallow",
                                                    "Halfway Exposed_Shallow",
                                                    "Red Beach_Shallow",
                                                    "Coral Gardens_Reef flat",
                                                    "Halfway_Reef flat",
                                                    "Coral Gardens_Shallow", 
                                                    "Coral Gardens_Deep",
                                                    "Halfway_Shallow",
                                                    "Halfway_Deep",
                                                    "First Point_Shallow",
                                                    "First Point_Deep",
                                                    "Eco_1_Shallow" ,
                                                    "Eco_1_Deep",
                                                    "Eco_2_Shallow" , 
                                                    "Eco_2_Deep",
                                                    "Blue Pools_BP"))

dat2$Morphology2 <- factor(dat2$Morphology2, levels = c("branched",
                                                        "unbranched",
                                                        "other"),
                                             labels = c("Branched",
                                                         "Unbranched",
                                                         "Other"))

str(dat2)

HKdatwI2021 <- dat2 %>% filter(Trip == "2021") %>% droplevels()
str(HKdatwI2021)

# Data without Blue Pools or First Point (Intermediate sites)
# Only exposed and sheltered 2021, surface layer.

HKdatwI2021NoBP <- HKdatwI2021 %>% dplyr::filter(sitedepth != "First Point_Shallow") %>%
  filter(sitedepth != "First Point_Deep") %>% 
  filter(sitedepth != "Blue Pools_BP") %>% droplevels()
levels(HKdatwI2021NoBP$sitedepth)
levels(HKdatwI2021NoBP$locexpdepth)
str(HKdatwI2021NoBP)

HKdatwI2021NoI <- HKdatwI2021 %>% dplyr::filter(sitedepth != "First Point_Shallow") %>%
  filter(sitedepth != "First Point_Deep") %>% droplevels()

levels(HKdatwI2021NoI$sitedepth)

names(HKdatwI2021NoBP)
str(HKdatwI2021NoBP)
levels(HKdatwI2021NoBP$Quadrat_angle2)

levels(HKdatwI2021NoI$Quadrat_angle2)

HKdatwI2021NoI %>% filter(Quadrat_angle2 != "NA") %>%
    ggplot(aes(x=locexpdepth, fill=Quadrat_angle2)) + 
    geom_bar(position = "fill")

# **Section 1** ----

# Habitat effects on rubble bed measurements (bed area, thickness)

# 1.1 Bed area ----

# Raw data (includes BP)

HKdatwI2021SUM_A <- HKdatwI2021 %>%
  filter(!is.na(Patch_area_m2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanarea = mean(Patch_area_m2, na.rm = TRUE),
            sdarea= sd(Patch_area_m2, na.rm = TRUE),
            searea = sdarea / sqrt(n()),
            ci95 = 1.96*searea) %>%
  mutate(lowerCI = meanarea - ci95,
         upperCI = meanarea + ci95) %>% as.data.frame()

view(HKdatwI2021SUM_A)

hist(log(HKdatwI2021NoBP$Patch_area_m2))

ggplot(HKdatwI2021SUM_A, aes(x = locexpdepth, y = meanarea)) +
  geom_point() + ylim(0,300)

bedsizeG <- glmmTMB(log(Patch_area_m2) ~ locexpdepth,
                       family = gaussian, data = HKdatwI2021NoBP)

car::Anova(bedsizeG)

plot(simulateResiduals(fittedModel = bedsizeG))

bedsizeGcomps <- emmeans(bedsizeG, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(bedsizeGcomps, "bedsizeGcomps_compsT.docx", "bedsizeGcomps_compsT2.docx")


# Plot

colsT <- c("Keppels (inshore)_Sheltered_Shallow" = "#6CA6CD",
           "Keppels (inshore)_Exposed_Shallow" = "#CD8C95",
           "Heron (offshore)_Sheltered_Reef flat" = "#AEE2F5",
           "Heron (offshore)_Sheltered_Shallow" = "#6CA6CD",
           "Heron (offshore)_Sheltered_Deep" = "#104E8B",
           "Heron (offshore)_Exposed_Shallow" = "#CD8C95",    
           "Heron (offshore)_Exposed_Deep" = "#A52A2A",
           "Dep. Area"  = "grey")



pd <- position_dodge(0.5)
legend_titleT <- "Habitat" 

# Raw data (plot BP raw data on plot)

HKdatwI2021SUM_A <- HKdatwI2021SUM_A %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_A$locexpdepth <- fct_recode(HKdatwI2021SUM_A$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

HKdatwI2021NoBP_grid <- with(HKdatwI2021NoBP, list(locexpdepth= levels(locexpdepth)))

bedsizeG.pred <- emmeans(bedsizeG, ~ locexpdepth, 
                          at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

(bedarea.plot <- ggplot() + 
    geom_point(data = bedsizeG.pred, aes(y=response, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bedsizeG.pred, aes(x = locexpdepth, ymin=lower.CL, ymax=upper.CL,
                                             colour = locexpdepth),
                  width=.2, position = pd) +
   ylim(0,2000) +
    labs(x="Habitat", y = "Mean rubble bed area (m2)") +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_A, aes(y= meanarea, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_A, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5.5 * 3.7


# 1.2 Bed thickness ----

hist(log(HKdatwI2021NoBP$Average_bed_depth_cm))
range(HKdatwI2021NoBP$Average_bed_depth_cm, na.rm=TRUE) # 2.3 to 51cm
mean(HKdatwI2021NoBP$Average_bed_depth_cm, na.rm=TRUE) # mean is 20cm
median(HKdatwI2021NoBP$Average_bed_depth_cm, na.rm=TRUE) # median is 17cm

bedthick <- glmmTMB(log(Average_bed_depth_cm) ~ locexpdepth + (1|sitedepth),
                       family = gaussian, data = HKdatwI2021NoBP)

car::Anova(bedthick)

plot(simulateResiduals(fittedModel = bedthick)) # decent

bedthickcomps <- emmeans(bedthick, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(bedthickcomps, "bedthickcomps_compsT.docx", "bedthickcomps_compsT2.docx")


# Raw data (includes BP)

HKdatwI2021SUM_Th <- HKdatwI2021 %>%
  filter(!is.na(Average_bed_depth_cm)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanthick = mean(Average_bed_depth_cm, na.rm = TRUE),
            sdthick= sd(Average_bed_depth_cm, na.rm = TRUE),
            sethick = sdthick / sqrt(n()),
            ci95 = 1.96*sethick) %>%
  mutate(lowerCI = meanthick - ci95,
         upperCI = meanthick + ci95) %>% as.data.frame()

# Plot

# Raw data (plot BP raw data on plot)

HKdatwI2021SUM_Th <- HKdatwI2021SUM_Th %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_Th$locexpdepth <- fct_recode(HKdatwI2021SUM_Th$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

bedthick.pred <- emmeans(bedthick, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

(bedthick.plot <- ggplot() + 
    geom_point(data = bedthick.pred, aes(y=response, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bedthick.pred, aes(x = locexpdepth, ymin=lower.CL, ymax=upper.CL,
                                            colour = locexpdepth),
                  width= .2, position = pd) +
    labs(x="Habitat", y = "Mean rubble bed thickness (cm)") +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_Th, aes(y= meanthick, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_Th, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                               colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5.5 * 3.7

# **Section 2** ----

# Habitat effects on rubble piece measurements

# 2.1 Void size -----

# Void size ~ habitat

hist(HKdatwI2021NoBP$Max_gap_size_to_right_cm) # very right skewed, do log  transform
HKdatwI2021NoBP$Max_gap_size_to_right_cm <- as.numeric(HKdatwI2021NoBP$Max_gap_size_to_right_cm)
range(HKdatwI2021NoBP$Max_gap_size_to_right_cm, na.rm=TRUE) # 0.1 to 10.2 cm

# Transform

hist(log(HKdatwI2021NoBP$Max_gap_size_to_right_cm)) # good

maxgap1 <- glmmTMB(log(Max_gap_size_to_right_cm2) ~ 
                     locexpdepth + (1|sitedepth/Quadrat2),
                   family = gaussian, data = HKdatwI2021NoBP) # don't want to investigate interaction due to difference in size between habitats, i.e., not representative sizes at all habitats.

car::Anova(maxgap1) # there is an interaction between habitat and rubble length on void size

# Model diagnostics

plot(simulateResiduals(fittedModel = maxgap1)) # looks good

maxgap1comps <- emmeans(maxgap1, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(maxgap1comps, "maxgap1comps_compsT.docx", "maxgap1comps_compsT2.docx")

# Plot

# Raw data (includes BP)

HKdatwI2021SUM_V <- HKdatwI2021 %>%
  filter(!is.na(Max_gap_size_to_right_cm2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meangap = mean(Max_gap_size_to_right_cm2, na.rm = TRUE),
            sdgap= sd(Max_gap_size_to_right_cm2, na.rm = TRUE),
            segap = sdgap / sqrt(n()),
            ci95 = 1.96*segap) %>%
  mutate(lowerCI = meangap - ci95,
         upperCI = meangap + ci95) %>% as.data.frame()


HKdatwI2021SUM_V <- HKdatwI2021SUM_V %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_V$locexpdepth <- fct_recode(HKdatwI2021SUM_V$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

maxgap1.pred <- emmeans(maxgap1, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

(maxgap1.plot <- ggplot() + 
    geom_point(data = maxgap1.pred, aes(y=response, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = maxgap1.pred, aes(x = locexpdepth, ymin=lower.CL, ymax=upper.CL,
                                            colour = locexpdepth),
                  width= .2, position = pd) +
    labs(x="Habitat", y = "Mean void size (cm)") +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_V, aes(y= meangap, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_V, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7


# 2.2 Rubble length ----

# Rubble length ~ habitat

HKdatwI2021NoBP$Widest_span_cm2 <- as.numeric(HKdatwI2021NoBP$Widest_span_cm2)

hist(HKdatwI2021NoBP$Widest_span_cm2) # right skewed - tail to right
range(HKdatwI2021NoBP$Widest_span_cm2, na.rm = TRUE) # 2.5cm up to 57cm
hist(log(HKdatwI2021NoBP$Widest_span_cm2)) # this looks good

rubsize1 <- glmmTMB(log(Widest_span_cm2) ~ 
                      locexpdepth + (1|sitedepth/Quadrat2),
                    family = gaussian, data = HKdatwI2021NoBP)
# random effect is quadrat within site - multiple rubble pieces per quadrat, & quadrat is nested in site

car::Anova(rubsize1)

rubsize1comps <- emmeans(rubsize1, pairwise ~ locexpdepth, type = "response") # 21 comparisons

save_tables_to_docx(rubsize1comps, "rubsize1comps_emmeans.docx", "rubsize1comps_contrasts.docx")


# Model diagnostics

plot(simulateResiduals(rubsize1)) # good

# Plot habitat

# Raw data (plot BP raw data on plot)

HKdatwI2021SUM_L <- HKdatwI2021 %>%
  filter(!is.na(Widest_span_cm2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanlength = mean(Widest_span_cm2, na.rm = TRUE),
            sdlength= sd(Widest_span_cm2, na.rm = TRUE),
            selength = sdlength / sqrt(n()),
            ci95 = 1.96*selength) %>%
  mutate(lowerCI = meanlength - ci95,
         upperCI = meanlength + ci95) %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_L$locexpdepth <- fct_recode(HKdatwI2021SUM_L$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")


rubsize1.pred <- emmeans(rubsize1, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

(rubsize1.plot <- ggplot() + 
    geom_point(data = rubsize1.pred, aes(y=response, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = rubsize1.pred, aes(x = locexpdepth, ymin=lower.CL, ymax=upper.CL,
                                            colour = locexpdepth),
                  width=.2, position = pd) +
    labs(x="Habitat", y = "Mean rubble length (cm)") +
    ylim(0,25) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_L, aes(y= meanlength, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_L, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                             colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7


# 2.3 Branchiness -----

# Branchiness ~ habitat

HKdatwI2021NoBP$Num_branches2 <- as.numeric(HKdatwI2021NoBP$Num_branches2)
HKdatwI2021NoBP$Num_branches3_0_for_other <- as.numeric(HKdatwI2021NoBP$Num_branches3_0_for_other)
str(HKdatwI2021NoBP)

brnum1 <- glmmTMB(Num_branches2 ~ 
                       locexpdepth + (1|sitedepth/Quadrat2), 
                     family = nbinom1(), data = HKdatwI2021NoBP)

brnum2 <- glmmTMB(Num_branches3_0_for_other ~ 
                    locexpdepth + (1|sitedepth/Quadrat2), 
                  family = nbinom1(), data = HKdatwI2021NoBP)

car::Anova(brnum1) # locexpdepth sig. - Chisq = 43.3; P < 0.001

car::Anova(brnum2) # locexpdepth sig. - Chisq = 39.7; P < 0.001


brnum1_comps <- emmeans(brnum1, pairwise ~ locexpdepth, type = "response") # 21 comparisons

save_tables_to_docx(brnum1_comps, "brnum1_comps_emmeans.docx", "brnum1_comps_contrasts.docx")


# So, generally, rubble in sheltered location at the same depth has more branches
# But that is not the case for Heron Sheltered Shallow vs Heron Exposed Shallow (P = 0.89,
# even though it's trending toward sheltered shallow having more branches - 2.6 vs 1.8)


emmeans(brnum2, pairwise ~ locexpdepth, type = "response") # 21 comparisons


# Model diagnostics

plot(simulateResiduals(brnum1)) # negative binomial is better than a poisson that was trialled


# Plot

# Raw data (plot BP raw data on plot)

HKdatwI2021SUM_Br <- HKdatwI2021 %>%
  filter(!is.na(Num_branches2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanbranch = mean(Num_branches2, na.rm = TRUE),
            sdbranch= sd(Num_branches2, na.rm = TRUE),
            sebranch = sdbranch / sqrt(n()),
            ci95 = 1.96*sebranch) %>%
  mutate(lowerCI = meanbranch - ci95,
         upperCI = meanbranch + ci95) %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_Br$locexpdepth <- fct_recode(HKdatwI2021SUM_Br$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")


# 2.4 Morphology -----

# Summary of % branched, unbranched, other:

HKdat1sum_BP <- HKdatwI2021NoI %>% group_by(locexpdepth, Morphology2) %>%
  summarise(N = sum(!is.na(Morphology2))) %>% # not including the counts of NAs
  group_by(locexpdepth) %>%
  mutate(Percentage = N/sum(N,na.rm = TRUE) * 100) %>% # again, not including the NAs
  filter(!is.na(Morphology2))

View(HKdat1sum_BP)
str(HKdatwI2021BP)

HKdat1sum_BPDF <- HKdat1sum_BP %>% as.data.frame()

write.csv(HKdat1sum_BPDF, "HKdat1sum_BPDF.csv")

# Plot

colsM <- c("Branched" = "#4F4F4F",
           "Unbranched" = "#999999",
           "Other" = "#E0E0E0")

legend_titleM <- "Morphology" 

(morphplot <- HKdatwI2021BP %>% filter(Morphology2 != "NA") %>%
    ggplot(aes(x=locexpdepth, fill=Morphology2)) + 
    geom_bar(position = "fill") +
    labs(x="Habitat", y = "Proportion of rubble pieces") +
    scale_fill_manual(legend_titleM, values = colsM, 
                      labels=c("Branched", "Unbranched","Other")) +
    scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="right") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5.5 * 3.7


# **Section 3** -----

# Habitat effects on binding & stability

# 3.1 Stability likelihood ----

# Look at raw data:

sumprop <- HKdatwI2021NoBP %>% 
  group_by(sitedepth) %>% summarise(stablemean = mean(Stable_0.1_v2, na.rm = TRUE))
# Halfway Shallow, 25% of pieces are not stable. 75% of pieces are stable.

sumprop <- HKdatwI2021NoBP %>% 
  group_by(sitedepth) %>% summarise(Boundmean = mean(Bound, na.rm = TRUE))
# Halfway Shallow, 25% of pieces are not stable. 75% of pieces are stable.

stabmod1 <- glmmTMB(Stable_0.1_v2 ~ 
                     locexpdepth +
                      (1|sitedepth/Quadrat2), family = "binomial", data = HKdatwI2021NoBP, # need to use different dataframe
                    na.action = "na.omit")  # Cannot look at morphology interactions - not enough representation of other and unbranched with size combos

plot(simulateResiduals(fittedModel = stabmod1)) # looks good

car::Anova(stabmod1) 

stabmod1comps <- emmeans(stabmod1, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(stabmod1comps, "stabmod1comps_emmeans.docx", "stabmod1comps_contrasts.docx")

# Raw data (plot BP raw data on plot)

HKdatwI2021SUM_Stab <- HKdatwI2021 %>%
  filter(!is.na(Stable_0.1_v2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanstab = mean(Stable_0.1_v2, na.rm = TRUE),
            sdstab= sd(Stable_0.1_v2, na.rm = TRUE),
            sestab = sdstab / sqrt(n()),
            ci95 = 1.96*sestab) %>%
  mutate(lowerCI = meanstab - ci95,
         upperCI = meanstab + ci95)

HKdatwI2021SUM_Stab <- HKdatwI2021SUM_Stab %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_Stab$locexpdepth <- fct_recode(HKdatwI2021SUM_Stab$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

stabmod1.pred <- emmeans(stabmod1, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

# Plot

(stabmod1.plot <- ggplot() + 
    geom_point(data = stabmod1.pred, aes(y=prob, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = stabmod1.pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                             colour = locexpdepth), width=.2, position = pd) +
    labs(x="Habitat", y = "Probability of stability on pick-up") +
    ylim(0,1) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_Stab, aes(y= meanstab, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_Stab, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                  colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7


# Stability types ----

# Stability types ~ habitat

# Summary of % branched, unbranched, other:

HKdatwI2021NoBP$Stable_description2 <- factor(HKdatwI2021NoBP$Stable_description2, 
                                              levels = c("unstable",
                                               "bound",
                                              "interlocked",
                                                "interlocked ",
                                              "buried",
                                               "heavy"),
                                              labels = c("Unstable",
                                              "Bound",
                                              "Interlocked",
                                               "Interlocked",
                                               "Buried",
                                                "Heavy"))

HKdatwI2021NoBPSUM <- HKdatwI2021NoBP %>% group_by(locexpdepth, Stable_description2) %>%
  summarise(N = sum(!is.na(Stable_description2))) %>% # not including the counts of NAs
  group_by(locexpdepth) %>%
  mutate(Percentage = N/sum(N,na.rm = TRUE) * 100) %>% # again, not including the NAs
  filter(!is.na(Stable_description2))

head(HKdatwI2021NoBPSUM)
levels(HKdatwI2021NoBP$Stable_description2)

# Model - Cumulative link

# relevel

HKdatwI2021NoBP_SD <- HKdatwI2021NoBP %>% 
  dplyr::filter(Stable_description2 != "Unstable") %>% droplevels()

HKdatwI2021NoBP_SD$Stable_description2 <- as.factor(HKdatwI2021NoBP_SD$Stable_description2)
levels(HKdatwI2021NoBP_SD$Stable_description2) # Unstable no longer included

HKdatwI2021NoBP_SD_SUM <- HKdatwI2021NoBP_SD %>% group_by(Stable_description2, locexpdepth) %>%
  summarise(N = sum(!is.na(Stable_description2))) 

HKdatwI2021NoBP_SD_SUM2 <- HKdatwI2021NoBP_SD %>% group_by(Stable_description2) %>%
  summarise(N = sum(!is.na(Stable_description2))) 

# Order the levels of Stable_description2 is bound, interlocked, buried and heavy, in order of usual stability

stabdes1 <- clmm(Stable_description2 ~ locexpdepth + (1|sitedepth), 
                    link = "logit", data = HKdatwI2021NoBP_SD, # can not include quadrat2 because of less reps now Unstable has been removed, some habitats only 1 rep of a certain level of stability (see HKdat1SDsum)
                    na.action = "na.omit") 
summary(stabdes1)

stabdes1comps <- emmeans(stabdes1, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(stabdes1comps, "stabdes1comp_compsT1s.docx", "stabdes1comps_compsT2.docx")

# Plot

HKdatwI2021NoI$Stable_description2 <- factor(HKdatwI2021NoI$Stable_description2, 
                                              levels = c("unstable",
                                                         "bound",
                                                         "interlocked",
                                                         "interlocked ",
                                                         "buried",
                                                         "heavy"),
                                              labels = c("Unstable",
                                                         "Bound",
                                                         "Interlocked",
                                                         "Interlocked",
                                                         "Buried",
                                                         "Heavy"))

HKdatwI2021BP_SD <- HKdatwI2021NoI %>% 
  filter(Stable_description2 != "Unstable") %>% droplevels()

head(HKdatwI2021BP_SD)

colsS <- c("Bound" = "#CCCCCC",
           "Interlocked" = "#949494",
           "Buried" = "#4F4F4F",
           "Heavy" = "#0A0A0A")

legend_titleM <- "Stability Mode" 

(stabilitydescplot <- HKdatwI2021BP_SD %>% filter(Stable_description2 != "NA") %>%
    ggplot(aes(x=locexpdepth, fill=Stable_description2)) + 
    geom_bar(position = "fill") +
    labs(x="Habitat", y = "Proportion of rubble pieces") +
    scale_fill_manual(legend_titleM, values = colsS, 
                      labels=c("Bound",
                               "Interlocked",
                               "Buried",
                               "Heavy")) +
    scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="right") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))


# 3.2 Binding likelihood -----

# With algae ----

bindmod1 <- glmmTMB(Bound ~ locexpdepth +
                      (1|sitedepth/Quadrat2), family = "binomial", data = HKdatwI2021NoBP, 
                    na.action = "na.omit")  # Cannot look at morphology interactions - not enough representation of other and unbranched with size combos

bindmodT <- glmmTMB(Bound ~ sitedepth +
                      (1|sitedepth/Quadrat2), family = "binomial", data = HKdatwI2021NoBP, 
                    na.action = "na.omit")  # Cannot look at morphology interactions - not enough representation of other and unbranched with size combos

car::Anova(bindmodT)

emmeans(bindmodT, pairwise ~ sitedepth, type = "response")

plot(simulateResiduals(fittedModel = bindmod1)) # good

car::Anova(bindmod1) 

bindmod1comps <- emmeans(bindmod1, pairwise  ~ locexpdepth, type = "response")

save_tables_to_docx(bindmod1comps, "bindmod1comps_emmeans.docx", "bindmod1comps_constrasts.docx")

#Raw data (plot BP raw data on plot)

HKdatwI2021SUM_Bind <- HKdatwI2021 %>%
  filter(!is.na(Bound)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanbound = mean(Bound, na.rm = TRUE),
            sdbound= sd(Bound, na.rm = TRUE),
            sebound = sdbound / sqrt(n()),
            ci95 = 1.96*sebound) %>%
  mutate(lowerCI = meanbound - ci95,
         upperCI = meanbound + ci95)

HKdatwI2021SUM_Bind <- HKdatwI2021SUM_Bind %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_Bind$locexpdepth <- fct_recode(HKdatwI2021SUM_Bind$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

bindmod1_pred <- emmeans(bindmod1, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

# Plot

(bindmod1.plot <- ggplot() + 
    geom_point(data = bindmod1_pred, aes(y=prob, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bindmod1_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                            colour = locexpdepth), width=.2, position = pd) +
    labs(x="Habitat", y = "Probability of binding") +
    ylim(0,1) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_Bind, aes(y= meanbound, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_Bind, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                  colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))

# Saved as 5 * 3.7


# Without algae ----

# Create dataset without algae binds

# Make dataset with Blue Pools so I can calculate the raw data, filter before running models

HKdat1_1 <- HKdatwI2021NoI %>% filter(!is.na(Bound))

# Then gather so in long form

HKdat1_long <- gather(HKdat1_1, key = "bindernumber", value = "type", 
                      Binder1.broad_cat, 
                      Binder2.broad_cat, 
                      Binder3.broad_cat, 
                      Binder4.broad_cat, 
                      Binder5.broad_cat, 
                      Binder6.broad_cat, 
                      Binder7.broad_cat, 
                      Binder8.broad_cat, 
                      Binder9.broad_cat,
                      Binder10.broad_cat, 
                      Binder11.broad_cat, 
                      Binder12.broad_cat, 
                      Binder13.broad_cat, 
                      Binder14.broad_cat, 
                      Binder15.broad_cat) # removing  na.rm=TRUE so that the zeroes are in there

# Now group by each piece, and if the type is NOT turf, macro, sum the number of rows for that piece.

View(HKdat1_long)

HKdat1_long$type <- as.factor(HKdat1_long$type)
levels(HKdat1_long$type)

# Make a new column where the turfs, cyanobacteria and macroalgae are NA as well.
# Then transform again to be 0 for NA and 1 for any other text and then sum

HKdat1_long2 <- HKdat1_long %>%
  mutate(typeNoAlgae = if_else(type %in% c(NA, "Cyanobacteria", "Lobophora", "Turf algae", "Other macroalgae"), NA_character_, type))

head(HKdat1_long2)

HKdat1_long3 <- HKdat1_long2 %>%
  mutate(typeNoAlgae = if_else(typeNoAlgae %in% c(NA), 0, 1)) # making cases 0 where they are NA

class(HKdat1_long3$typeNoAlgae)
class(HKdat1_long3$Unique_piece)

HKdat1_long3SUM <- HKdat1_long3 %>% group_by(Unique_piece,
                                             Trip, Location, Exposure, Depth, Site,                                      
                                             Quadrat2, Patch_area_m2, Average_bed_depth_cm,Piece,                                    
                                             Number_contact_points2, Max_gap_size_to_right_cm2, 
                                             Stable_0.1_v2, Stable_description2,                          
                                             Widest_span_cm2, Morphology2, Num_branches2,Bound, Number_binds2,
                                             Bound_to_how_many,exposuredepth, locexpdepth, sitedepth) %>%
  summarise(totalbinds_NoA = sum(as.numeric(typeNoAlgae), na.rm = TRUE),
            boundo_NoA = max(as.numeric(typeNoAlgae), na.rm = TRUE))

str(HKdat1_long3SUM)
levels(HKdat1_long3SUM$sitedepth)
# Was 780 observations without the Bound = NAs removed. Now it is 762 (without Blue Pools)
# Some rubble pieces not assessed in quadrats due to time running out underwater

HKdat1_long3SUM_NoBP <- HKdat1_long3SUM %>% filter(sitedepth != "Blue Pools_BP") %>% droplevels()
levels(HKdat1_long3SUM_NoBP$sitedepth)

bindmod1NA <- glmmTMB(boundo_NoA ~ locexpdepth +
                        (1|sitedepth/Quadrat2), family = "binomial", data = HKdat1_long3SUM_NoBP, 
                      na.action = "na.omit")

car::Anova(bindmod1NA)

plot(simulateResiduals(fittedModel = bindmod1NA)) # good

bindmod1NAcomps <- emmeans(bindmod1NA, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(bindmod1NAcomps, "bindmod1NAcomps_emmeans.docx", "bindmod1NAcomps_constrasts.docx")

# Plot 

#Raw data (plot BP raw data on plot)

HKdatwI2021SUM_BoundNoA <- HKdat1_long3SUM %>%
  filter(!is.na(boundo_NoA)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanboundNoA = mean(boundo_NoA, na.rm = TRUE),
            sdboundNoA= sd(boundo_NoA, na.rm = TRUE),
            seboundNoA = sdboundNoA / sqrt(n()),
            ci95 = 1.96*seboundNoA) %>%
  mutate(lowerCI = meanboundNoA - ci95,
         upperCI = meanboundNoA + ci95)

HKdatwI2021SUM_BoundNoA <- HKdatwI2021SUM_BoundNoA %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_BoundNoA$locexpdepth <- fct_recode(HKdatwI2021SUM_BoundNoA$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

bindmod1NA_pred <- emmeans(bindmod1NA, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

# Plot

(bindmod1NA_plot <- ggplot() + 
    geom_point(data = bindmod1NA_pred, aes(y=prob, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bindmod1NA_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                            colour = locexpdepth), width=.2, position = pd) +
    labs(x="Habitat", y = "Probability of binding (No algae)") +
    ylim(0,1) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_BoundNoA, aes(y= meanboundNoA, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_BoundNoA, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                       colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7

# Plot with both

(bindmod1Both_plot <- ggplot() + 
   geom_point(data = bindmod1NA_pred, aes(y=prob, x=locexpdepth, colour = locexpdepth), 
              position=pd, size =3) +
   geom_errorbar(data = bindmod1NA_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                             colour = locexpdepth), width=.2, position = pd) +
   geom_point(data = bindmod1_pred, aes(y=prob, x=locexpdepth, colour = locexpdepth), 
              position=pd, size =3) +
   geom_errorbar(data = bindmod1_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                           colour = locexpdepth), width=.2, position = pd) +
   labs(x="Habitat", y = "Probability of binding") +
   ylim(0,1) +
   scale_colour_manual(legend_titleT, values=colsT, 
                       labels=c("In Shelt Shal","In Exp Shal",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
   scale_fill_manual(legend_titleT, values = colsT, 
                     labels=c("In Shelt Shal","In Exp Shal",
                              "Off Shelt Shal", "Off Shelt Deep", 
                              "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
   scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                    labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
   geom_point(data = HKdatwI2021SUM_BoundNoA, aes(y= meanboundNoA, x = locexpdepth, colour = locexpdepth), 
              position=pd, size =3) + 
   geom_errorbar(data = HKdatwI2021SUM_BoundNoA, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                     colour = locexpdepth), 
                 width=.2, position = pd) +
   geom_point(data = HKdatwI2021SUM_Bind, aes(y= meanbound, x = locexpdepth, colour = locexpdepth), 
              position=pd, size =3) + 
   geom_errorbar(data = HKdatwI2021SUM_Bind, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                 colour = locexpdepth), 
                 width=.2, position = pd) +
   theme_minimal_hgrid(line_size = 0.2) +
   theme(text = element_text(size=10, face="bold")) + 
   theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
   theme(panel.background =element_rect(colour = "black", size=1)) + 
   theme(axis.ticks.length=unit(.2,"cm")) +
   theme(legend.position="none") +
   theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7

# 3.3 Binding density ----

# With algae -----

hist(HKdatwI2021NoBP$Number_binds2) # poisson or neg. binomial

bindtot1 <- glmmTMB(Number_binds2 ~ locexpdepth + 
                           (1|sitedepth/Quadrat2), family = poisson, data = HKdatwI2021NoBP, 
                         na.action = "na.omit")

car::Anova(bindtot1)

plot(simulateResiduals(fittedModel = bindtot1)) # good

bindtot1comps <- emmeans(bindtot1, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(bindtot1comps, "bindtot1comps_emmeans.docx", "bindtot1comps_contrasts.docx")

#Raw data (plot BP raw data on plot)

HKdatwI2021SUM_TotBinds <- HKdatwI2021 %>%
  filter(!is.na(Number_binds2)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meanbinds = mean(Number_binds2, na.rm = TRUE),
            sdbinds= sd(Number_binds2, na.rm = TRUE),
            sebinds = sdbinds / sqrt(n()),
            ci95 = 1.96*sebinds) %>%
  mutate(lowerCI = meanbinds - ci95,
         upperCI = meanbinds + ci95)

HKdatwI2021SUM_TotBindsv <- HKdatwI2021SUM_TotBinds %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_TotBindsv$locexpdepth <- fct_recode(HKdatwI2021SUM_TotBindsv$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

bindtot1_pred <- emmeans(bindtot1, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

# Plot

(bindtot1_plot <- ggplot() + 
    geom_point(data = bindtot1_pred, aes(y=rate, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bindtot1_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                            colour = locexpdepth), width=.2, position = pd) +
    labs(x="Habitat", y = "Number of binds") +
    ylim(0,3.3) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_TotBindsv, aes(y= meanbinds, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_TotBindsv, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                  colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7


# Without algae -----

names(HKdat1_long3SUM)


bindtot1NA <- glmmTMB(totalbinds_NoA ~ locexpdepth +
                        (1|sitedepth/Quadrat2), family = poisson, data = HKdat1_long3SUM_NoBP, 
                      na.action = "na.omit")

car::Anova(bindtot1NA)

plot(simulateResiduals(fittedModel = bindtot1NA)) # good

bindtot1NAcomps <- emmeans(bindtot1NA, pairwise ~ locexpdepth, type = "response")

save_tables_to_docx(bindtot1NAcomps, "bindtot1NAcomps_emmeans.docx", "bindtot1NAcomps_contrasts.docx")

# Plot

#Raw data (plot BP raw data on plot)

HKdatwI2021SUM_TotBindsNoA <- HKdat1_long3SUM %>%
  filter(!is.na(totalbinds_NoA)) %>% # or filter here to remove the NAs first
  group_by(locexpdepth) %>%
  summarise(meantotbindsNoA = mean(totalbinds_NoA, na.rm = TRUE),
            sdtotbindsNoA= sd(totalbinds_NoA, na.rm = TRUE),
            setotbindsNoA = sdtotbindsNoA / sqrt(n()),
            ci95 = 1.96*setotbindsNoA) %>%
  mutate(lowerCI = meantotbindsNoA - ci95,
         upperCI = meantotbindsNoA + ci95)

HKdatwI2021SUM_TotBindsNoA <- HKdatwI2021SUM_TotBindsNoA %>% as.data.frame() %>%
  filter(locexpdepth == "Heron (offshore)_Intermediate_BP") %>% droplevels()

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_TotBindsNoA$locexpdepth <- fct_recode(HKdatwI2021SUM_TotBindsNoA$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")


bindtot1NA_pred <- emmeans(bindtot1NA, ~ locexpdepth, 
                         at=HKdatwI2021NoBP_grid, type = 'response') %>% as.data.frame

# Plot

(bindtot1NA_plot <- ggplot() + 
    geom_point(data = bindtot1NA_pred, aes(y=rate, x=locexpdepth, colour = locexpdepth), 
               position=pd, size =3) +
    geom_errorbar(data = bindtot1NA_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                            colour = locexpdepth), width=.2, position = pd) +
    labs(x="Habitat", y = "Number of binds (No algae)") +
    ylim(0,3.3) +
    scale_colour_manual(legend_titleT, values=colsT, 
                        labels=c("In Shelt Shal","In Exp Shal",
                                 "Off Shelt Shal", "Off Shelt Deep", 
                                 "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_fill_manual(legend_titleT, values = colsT, 
                      labels=c("In Shelt Shal","In Exp Shal",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
    scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                     labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    geom_point(data = HKdatwI2021SUM_TotBindsNoA, aes(y= meantotbindsNoA, x = locexpdepth, colour = locexpdepth), 
               position=pd, size =3) + 
    geom_errorbar(data = HKdatwI2021SUM_TotBindsNoA, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                       colour = locexpdepth), 
                  width=.2, position = pd) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7

# Total Plot ----

(bindtotandbound1_plot <- ggplot() + 
   geom_point(data = bindtot1_pred, aes(y=rate, x=locexpdepth, colour = locexpdepth), 
              position=pd, size =3) +
   geom_errorbar(data = bindtot1_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                           colour = locexpdepth), width=.2, position = pd) +
   labs(x="Habitat", y = "Number of binds") +
   geom_point(data = bindtot1NA_pred, aes(y=rate, x=locexpdepth, colour = locexpdepth), 
              position=pd, size =3) +
   geom_errorbar(data = bindtot1NA_pred, aes(x = locexpdepth, ymin=asymp.LCL, ymax=asymp.UCL,
                                             colour = locexpdepth), width=.2, position = pd) +
   ylim(0,3.3) +
   geom_point(data = HKdatwI2021SUM_TotBindsNoA, aes(y= meantotbindsNoA, x = locexpdepth, colour = locexpdepth), 
              position=pd, size =3) + 
   geom_errorbar(data = HKdatwI2021SUM_TotBindsNoA, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                        colour = locexpdepth), 
                 width=.2, position = pd) +
   scale_colour_manual(legend_titleT, values=colsT, 
                       labels=c("In Shelt Shal","In Exp Shal",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
   scale_fill_manual(legend_titleT, values = colsT, 
                     labels=c("In Shelt Shal","In Exp Shal",
                              "Off Shelt Shal", "Off Shelt Deep", 
                              "Off Exp Shal", "Off Exp Deep", "Off Shelt RF")) +
   scale_x_discrete(limits = c(levels(HKdatwI2021NoBP$locexpdepth), "Dep. Area"),
                    labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                               "Off Shelt Shal", "Off Shelt Deep", 
                               "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
   geom_point(data = HKdatwI2021SUM_TotBindsv, aes(y= meanbinds, x = locexpdepth, colour = locexpdepth), 
              position=pd, size =3) + 
   geom_errorbar(data = HKdatwI2021SUM_TotBindsv, aes(x = locexpdepth, ymin=lowerCI, ymax=upperCI, 
                                                      colour = locexpdepth), 
                 width=.2, position = pd) +
   theme_minimal_hgrid(line_size = 0.2) +
   theme(text = element_text(size=10, face="bold")) + 
   theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
   theme(panel.background =element_rect(colour = "black", size=1)) + 
   theme(axis.ticks.length=unit(.2,"cm")) +
   theme(legend.position="none") +
   theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5 * 3.7

# **Section 4** -----

# Rubble bed & piece effects on binding & stability
# random effect of habitat/site/quadrat

# 4.1 Stability likelihood ----

# Stability likelihood ~ bed and piece characteristics
# Interactions will not run.

stabmod1R <- glmmTMB(Stable_0.1_v2 ~ 
                         Widest_span_cm2 + 
                         Morphology2 + 
                         Num_branches2 +
                         Max_gap_size_to_right_cm2 +
                         Average_bed_depth_cm +
                         scale(Patch_area_m2) +
                         (1|locexpdepth/sitedepth/Quadrat2), family = "binomial", 
                       data = HKdatwI2021NoBP,
                     na.action = "na.omit") 

car::Anova(stabmod1R) # only widest span is important
plot(simulateResiduals(fittedModel = stabmod1R))

# Remove non-sig. 

stabmod2R <- update(stabmod1R, ~ . -  scale(Patch_area_m2))
car::Anova(stabmod2R) 
AICc(stabmod1R, stabmod2R)

stabmod3R <- update(stabmod2R, ~ . -  Average_bed_depth_cm)
car::Anova(stabmod3R) 
AICc(stabmod2R, stabmod3R)# AIC goes up, leave in

stabmod4R <- update(stabmod2R, ~ . -  Max_gap_size_to_right_cm2)
car::Anova(stabmod4R) 
AICc(stabmod2R, stabmod4R)

stabmod5R <- update(stabmod4R, ~ . -  Morphology2)
car::Anova(stabmod5R) 
AICc(stabmod4R, stabmod5R)

stabmod6R <- update(stabmod5R, ~ . -  Num_branches2)
car::Anova(stabmod6R) 
AICc(stabmod5R, stabmod6R) # don't remove Num branches

stabmod5Rcsv <- car::Anova(stabmod5R) %>% as.data.frame()

write.csv(stabmod5Rcsv, "stabmod5R.csv")

summary(stabmod5R)

# Plot Rubble length

length.grid <- with(HKdatwI2021NoBP, 
                         list(Widest_span_cm2 = seq(min(Widest_span_cm2, na.rm = TRUE), 
                                                    max(Widest_span_cm2, na.rm = TRUE),len = 100)))

stabmod5R.pred <- emmeans(stabmod5R, ~ Widest_span_cm2, 
                          at=length.grid, type = 'response') %>% as.data.frame()
head(stabmod5R.pred)

(stabmod5R.plot <- ggplot(stabmod5R.pred, aes(y=prob, x=Widest_span_cm2)) +
    geom_point(aes(x = Widest_span_cm2, y = Stable_0.1_v2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble length (cm)", y = "Probability of stability") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# emmeans

lengths <- with(HKdatwI2021NoBP, 
                    list(Widest_span_cm2 = c(3, 10, 20)))

stabmod5Rcomps <- emmeans(stabmod5R, pairwise ~ Widest_span_cm2, at = lengths, type = "response")


# Plot Average bed depth

beddepth.grid <- with(HKdatwI2021NoBP, 
                    list(Average_bed_depth_cm = seq(min(Average_bed_depth_cm, na.rm = TRUE), 
                                               max(Average_bed_depth_cm, na.rm = TRUE),len = 100)))

stabmod5R.pred2 <- emmeans(stabmod5R, ~ Average_bed_depth_cm, 
                           at=beddepth.grid, type = 'response') %>% as.data.frame()
head(stabmod5R.pred2)

(stabmod5R.plot <- ggplot(stabmod5R.pred2, aes(y=prob, x=Average_bed_depth_cm)) +
    geom_point(aes(x = Average_bed_depth_cm, y = Stable_0.1_v2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble bed thickness (cm)", y = "Probability of stability") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# emmeans

depths <- with(HKdatwI2021NoBP, 
                list(Average_bed_depth_cm = c(2, 20, 35)))

stabmod5Rcomps2 <- emmeans(stabmod5R, pairwise ~ Average_bed_depth_cm, at = depths, type = "response")


# 4.2 Binding likelihood ----

# With algae ----

bindmod1R <- glmmTMB(Bound ~ 
                       Widest_span_cm2 + 
                       Morphology2 + 
                       Num_branches2 +
                       Max_gap_size_to_right_cm2 +
                       Average_bed_depth_cm +
                       scale(Patch_area_m2) +
                       (1|locexpdepth/sitedepth/Quadrat2), family = "binomial", 
                     data = HKdatwI2021NoBP,
                     na.action = "na.omit") 

plot(simulateResiduals(fittedModel = bindmod1R)) # good

car::Anova(bindmod1R)

# Remove non-sig. 

bindmod2R <- update(bindmod1R, ~ . -  scale(Patch_area_m2))
car::Anova(bindmod2R) 
AICc(bindmod1R, bindmod2R)

bindmod3R <- update(bindmod2R, ~ . -  Average_bed_depth_cm)
car::Anova(bindmod3R) 
AICc(bindmod2R, bindmod3R)# AIC goes up, leave in

bindmod4R <- update(bindmod2R, ~ . -  Max_gap_size_to_right_cm2)
car::Anova(bindmod4R) 
AICc(bindmod2R, bindmod4R)

bindmod5R <- update(bindmod4R, ~ . -  Morphology2)
car::Anova(bindmod5R) 
AICc(bindmod4R, bindmod5R)

bindmod6R <- update(bindmod5R, ~ . -  Num_branches2)
car::Anova(bindmod6R) 
AICc(bindmod5R, bindmod6R) # don't remove Num branches

car::Anova(bindmod5R)
summary(bindmod5R)

# Plot Rubble Length

bindmod5R_pred <- emmeans(bindmod5R, ~ Widest_span_cm2, 
                            at=length.grid, type = 'response') %>% as.data.frame()
head(bindmod5R_pred)

(bindmod5R_plot <- ggplot(bindmod5R_pred, aes(y=prob, x=Widest_span_cm2)) +
    geom_point(aes(x = Widest_span_cm2, y = Bound), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble widest span (cm)", y = "Probability of rubble being bound") +
    ylim(0,1) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))

# Plot Average bed depth

bindmod5R_pred2 <- emmeans(bindmod5R, ~ Average_bed_depth_cm, 
                             at=beddepth.grid, type = 'response') %>% as.data.frame

(bindmod5R_plot2 <- ggplot(bindmod5R_pred2, aes(y=prob, x=Average_bed_depth_cm)) +
    geom_point(aes(x = Average_bed_depth_cm, y = Bound), data = HKdatwI2021BP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble bed thickness (cm)", y = "Probability of rubble being bound") +
    ylim(0,1) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))

# Without algae ----

HKdat1_long3SUM_NoBP$Num_branches2 <- as.numeric(HKdat1_long3SUM_NoBP$Num_branches2)

bindprob1NA_R <- glmmTMB(boundo_NoA ~ Widest_span_cm2 + 
                        Morphology2 + 
                        Num_branches2 +
                        Max_gap_size_to_right_cm2 +
                        Average_bed_depth_cm +
                        scale(Patch_area_m2) +
                        (1|locexpdepth/sitedepth/Quadrat2), family = "binomial", 
                      data = HKdat1_long3SUM_NoBP,
                      na.action = "na.omit") 

car::Anova(bindprob1NA_R)

plot(simulateResiduals(fittedModel = bindprob1NA_R)) # good

bindprob2NA_R <- update(bindprob1NA_R, ~ . -  scale(Patch_area_m2))
car::Anova(bindprob2NA_R) 
AICc(bindprob1NA_R, bindprob2NA_R)

bindprob3NA_R <- update(bindprob2NA_R, ~ . -  Average_bed_depth_cm)
car::Anova(bindprob3NA_R) 
AICc(bindprob2NA_R, bindprob3NA_R)# AIC goes up, leave in

bindprob3NA_R <- update(bindprob2NA_R, ~ . -  Max_gap_size_to_right_cm2)
car::Anova(bindprob3NA_R) 
AICc(bindprob2NA_R, bindprob3NA_R)

bindprob4NA_R <- update(bindprob3NA_R, ~ . -  Morphology2)
car::Anova(bindprob4NA_R) 
AICc(bindprob3NA_R, bindprob4NA_R)

bindprob5NA_R <- update(bindprob4NA_R, ~ . -  Num_branches2)
car::Anova(bindprob5NA_R) 
AICc(bindprob4NA_R, bindprob5NA_R) # AIC goes up, don't remove Num branches

car::Anova(bindprob4NA_R) # trend toward average bed depth having effect.
# Widest span still has an affect for non-algae binds only.

summary(bindprob5NA_R)

# Plot 

bindprob4NA_R_pred <- emmeans(bindprob4NA_R, ~ Widest_span_cm2, 
                           at=length.grid2, type = 'response') %>% as.data.frame

# Plot

(bindprob4NA_R_plot <- ggplot(bindprob4NA_R_pred, aes(y=prob, x=Widest_span_cm2)) +
    geom_point(aes(x = Widest_span_cm2, y = Bound), data = HKdat1_long3SUM_NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble length (cm)", y = "Probability of binding (no algae)") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# Plot 

bindprob4NA_R_pred2 <- emmeans(bindprob4NA_R, ~ Average_bed_depth_cm, 
                              at=beddepth.grid2, type = 'response') %>% as.data.frame

# Plot

(bindprob4NA_R_plot2 <- ggplot(bindprob4NA_R_pred2, aes(y=prob, x=Average_bed_depth_cm)) +
    geom_point(aes(x = Average_bed_depth_cm, y = Bound), data = HKdat1_long3SUM_NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    labs(x="Rubble bed thickness (cm)", y = "Probability of binding (no algae)") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# Plots with algae AND not algae on same plot

(bindmod_R_plot_BOTH <- ggplot() + 
  geom_point(data = bindmod5R_pred, aes(y=prob, x=Widest_span_cm2), 
             position=pd, size =0.5, alpha= 1) +
            geom_ribbon(data = bindmod5R_pred, aes(x = Widest_span_cm2, ymin=asymp.LCL, 
                          ymax=asymp.UCL), alpha =0.3) +
            geom_point(data = bindprob4NA_R_pred, aes(y=prob, x=Widest_span_cm2), 
                       colour = "#CBB6F5", position=pd, size =0.5, alpha= 2) +
            geom_ribbon(data = bindprob4NA_R_pred, aes(x = Widest_span_cm2, ymin=asymp.LCL, 
                        ymax=asymp.UCL), fill = "#CBB6F5", 
                        alpha =0.3) +
    geom_point(aes(x = Widest_span_cm2, y = Bound), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_point(aes(x = Widest_span_cm2, y = boundo_NoA), data = HKdat1_long3SUM_NoBP, colour ="#CBB6F5", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
             labs(x="Rubble length (cm)", y = "Probability of binding") +
             ylim(0,1) +
            theme_minimal_hgrid(line_size = 0.2) +
            theme_minimal_hgrid(line_size = 0.2) +
            theme(text = element_text(size=10, face="bold")) + 
           theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
            theme(panel.background =element_rect(colour = "black", size=1)) + 
            theme(axis.ticks.length=unit(.2,"cm")) +
            theme(legend.position="none") +
            theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# Plots with algae AND not algae on same plot - Bed thickness

(bindmod_R_plot_BOTH2 <- ggplot() + 
    geom_point(data = bindmod5R_pred2, aes(y=prob, x=Average_bed_depth_cm), 
               position=pd, size =0.5, alpha= 1) +
    geom_ribbon(data = bindmod5R_pred2, aes(x = Average_bed_depth_cm, ymin=asymp.LCL, 
                                           ymax=asymp.UCL), alpha =0.3) +
    geom_point(data = bindprob4NA_R_pred2, aes(y=prob, x=Average_bed_depth_cm), 
               colour = "#CBB6F5", position=pd, size =0.5, alpha= 2) +
    geom_ribbon(data = bindprob4NA_R_pred2, aes(x = Average_bed_depth_cm, ymin=asymp.LCL, 
                                               ymax=asymp.UCL), fill = "#CBB6F5", 
                alpha =0.3) +
    geom_point(aes(x = Average_bed_depth_cm, y = Bound), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_point(aes(x = Average_bed_depth_cm, y = boundo_NoA), data = HKdat1_long3SUM_NoBP, colour ="#CBB6F5", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    labs(x="Rubble bed thickness (cm)", y = "Probability of binding") +
    ylim(0,1) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# 4.3 Binding density ----

# With algae -----

bindtot1R <- glmmTMB(Number_binds2 ~ Widest_span_cm2 + 
                      Morphology2 + 
                      Num_branches2 +
                      Max_gap_size_to_right_cm2 +
                      Average_bed_depth_cm +
                      scale(Patch_area_m2) +
                      (1|locexpdepth/sitedepth/Quadrat2), family = poisson, 
                    data = HKdatwI2021NoBP,
                    na.action = "na.omit") #nbinom2 didn't run,use poisson

car::Anova(bindtot1R) # only widest span important to number of binds

bindtot2R <- update(bindtot1R, ~ . - scale(Patch_area_m2))
AICc(bindtot1R, bindtot2R) # AIC decreases
car::Anova(bindtot1R)

bindtot3R <- update(bindtot2R, ~ . - Max_gap_size_to_right_cm2)
AICc(bindtot2R, bindtot3R) # AIC decreases
car::Anova(bindtot3R)

bindtot4R <- update(bindtot3R, ~ . - Num_branches2)
AICc(bindtot3R, bindtot4R) # AIC increases, leave in
car::Anova(bindtot4R)

bindtot4R <- update(bindtot3R, ~ . - Morphology2)
AICc(bindtot3R, bindtot4R) # AIC decreases
car::Anova(bindtot4R)

summary(bindtot4R)

plot(simulateResiduals(fittedModel = bindtot4R)) # decent

# Plot widest span

bindtot4R.pred <- emmeans(bindtot4R, ~ Widest_span_cm2, 
                               at=length.grid, type = 'response') %>% as.data.frame()
head(bindtot4R.pred)

(bindtot4R.plot <- ggplot(bindtot4R.pred, aes(y=rate, x=Widest_span_cm2)) +
    geom_point(aes(x = Widest_span_cm2, y = Number_binds2), data = HKdatwI2021BP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    ylim(0,8) +
    labs(x="Rubble length (cm)", y = "Number of binds per piece") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))

# Bed depth trending toward having impact on the number of binds

# Without algae -----

bindtot1NA_R <- glmmTMB(totalbinds_NoA ~ Widest_span_cm2 + 
                          Morphology2 + 
                          Num_branches2 +
                          Max_gap_size_to_right_cm2 +
                          Average_bed_depth_cm +
                          scale(Patch_area_m2) +
                          (1|locexpdepth/sitedepth/Quadrat2), family = poisson, 
                          data = HKdat1_long3SUM_NoBP, 
                          na.action = "na.omit")

car::Anova(bindtot1NA_R)


bindtot2NA_R <- update(bindtot1NA_R, ~ . - scale(Patch_area_m2))
AICc(bindtot1NA_R, bindtot2NA_R) # AIC decreases
car::Anova(bindtot2NA_R)

bindtot3NA_R <- update(bindtot2NA_R, ~ . - Max_gap_size_to_right_cm2)
AICc(bindtot2NA_R, bindtot3NA_R) # AIC decreases
car::Anova(bindtot3NA_R)

bindtot4NA_R <- update(bindtot3NA_R, ~ . - Morphology2)
AICc(bindtot3NA_R, bindtot4NA_R) # AIC decreases
car::Anova(bindtot4NA_R)
summary(bindtot4NA_R)

# Num branches, widest span & average bed depth ALL dig.
# Whereas when algae was included,only widest span was sig. with a trend toward bed depth

plot(simulateResiduals(fittedModel = bindtot4NA_R)) # decent

# Plot Rubble length (algae and no algae)

# Change to calling it Dep. Area so I can plot it on the same plot below.
HKdatwI2021SUM_TotBindsNoA$locexpdepth <- fct_recode(HKdatwI2021SUM_TotBindsNoA$locexpdepth, "Dep. Area" = "Heron (offshore)_Intermediate_BP")

length.grid2 <- with(HKdat1_long3SUM_NoBP, 
                    list(Widest_span_cm2 = seq(min(Widest_span_cm2, na.rm = TRUE), 
                                             max(Widest_span_cm2, na.rm = TRUE),len = 100)))

bindtot4NA_R_predL <- emmeans(bindtot4NA_R, ~ Widest_span_cm2, 
                              at=length.grid2, type = 'response') %>% as.data.frame


(bindtot_R_plot_BOTH_Length <- ggplot() + 
    geom_point(data = bindtot4R.pred, aes(y=rate, x=Widest_span_cm2), 
               position=pd, size =0.5, alpha= 1) +
    geom_ribbon(data = bindtot4R.pred, aes(x = Widest_span_cm2, ymin=asymp.LCL, 
                                           ymax=asymp.UCL), alpha =0.3) +
    geom_point(data = bindtot4NA_R_predL, aes(y=rate, x=Widest_span_cm2), 
               colour = "#CBB6F5", position=pd, size =0.5, alpha= 2) +
    geom_ribbon(data = bindtot4NA_R_predL, aes(x = Widest_span_cm2, ymin=asymp.LCL, 
                                               ymax=asymp.UCL), fill = "#CBB6F5", 
                alpha =0.3) +
    ylim(0,6) +
    geom_point(aes(x = Widest_span_cm2, y = Number_binds2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_point(aes(x = Widest_span_cm2, y = totalbinds_NoA), data = HKdat1_long3SUM_NoBP, colour ="#CBB6F5", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    labs(x="Rubble length (cm)", y = "Number of binds") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07

# Plot Average bed depth (Algae and no algae)

beddepth.grid2 <- with(HKdat1_long3SUM_NoBP, 
                    list(Average_bed_depth_cm = seq(min(Average_bed_depth_cm, na.rm = TRUE), 
                                             max(Average_bed_depth_cm, na.rm = TRUE),len = 100)))

bindtot4NA_R_predD <- emmeans(bindtot4NA_R, ~ Average_bed_depth_cm, 
                              at=beddepth.grid2, type = 'response') %>% as.data.frame

bindtot4R.pred_Depth <- emmeans(bindtot4R, ~ Average_bed_depth_cm, 
                          at=beddepth.grid, type = 'response') %>% as.data.frame()

(bindtot_R_plot_BOTH_Depth <- ggplot() + 
    geom_point(data = bindtot4R.pred_Depth, aes(y=rate, x=Average_bed_depth_cm), 
               position=pd, size =0.5, alpha= 1) +
    geom_ribbon(data = bindtot4R.pred_Depth, aes(x = Average_bed_depth_cm, ymin=asymp.LCL, 
                                           ymax=asymp.UCL), alpha =0.3) +
    geom_point(data = bindtot4NA_R_predD, aes(y=rate, x=Average_bed_depth_cm), 
               colour = "#CBB6F5", position=pd, size =0.5, alpha= 2) +
    geom_ribbon(data = bindtot4NA_R_predD, aes(x = Average_bed_depth_cm, ymin=asymp.LCL, 
                                               ymax=asymp.UCL), fill = "#CBB6F5", 
                alpha =0.3) +
    ylim(0,4)+
    geom_point(aes(x = Average_bed_depth_cm, y = Number_binds2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_point(aes(x = Average_bed_depth_cm, y = totalbinds_NoA), data = HKdat1_long3SUM_NoBP, colour ="#CBB6F5", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    labs(x="Rubble bed thickness (cm)", y = "Number of binds") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07


# Plot Number of branches (Algae and no algae)

branch.grid <- with(HKdatwI2021NoBP, 
                    list(Num_branches2 = seq(min(Num_branches2, na.rm = TRUE), 
                                               max(Num_branches2, na.rm = TRUE),len = 100)))

bindtot4R.pred <- emmeans(bindtot4R, ~ Num_branches2, 
                          at=branch.grid, type = 'response') %>% as.data.frame()

branch.grid2 <- with(HKdat1_long3SUM_NoBP, 
                    list(Num_branches2 = seq(min(Num_branches2, na.rm = TRUE), 
                                             max(Num_branches2, na.rm = TRUE),len = 100)))

bindtot4NA_R_predBr <- emmeans(bindtot4NA_R, ~ Num_branches2, 
                               at=branch.grid2, type = 'response') %>% as.data.frame

(bindtot_R_plot_BOTH_Br <- ggplot() + 
    geom_point(data = bindtot4R.pred, aes(y=rate, x=Num_branches2), 
               position=pd, size =0.5, alpha= 1) +
    geom_ribbon(data = bindtot4R.pred, aes(x = Num_branches2, ymin=asymp.LCL, 
                                                 ymax=asymp.UCL), alpha =0.3) +
    geom_point(data = bindtot4NA_R_predBr, aes(y=rate, x=Num_branches2), 
               colour = "#CBB6F5", position=pd, size =0.5, alpha= 2) +
    geom_ribbon(data = bindtot4NA_R_predBr, aes(x = Num_branches2, ymin=asymp.LCL, 
                                               ymax=asymp.UCL), fill = "#CBB6F5", 
                alpha =0.3) +
    ylim(0,2)+
    xlim(0,9.5) +
    geom_point(aes(x = Num_branches2, y = Number_binds2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_point(aes(x = Num_branches2, y = totalbinds_NoA), data = HKdat1_long3SUM_NoBP, colour ="#CBB6F5", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    labs(x="Branchiness", y = "Number of binds") +
    theme_minimal_hgrid(line_size = 0.2) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# Save as 5 * 3.07


# **Section 5** ----

# Predicting things using rubble length


# 5.1 Bed thickness ----

ggplot(HKdatwI2021NoBP, aes(y = Average_bed_depth_cm, x = Widest_span_cm2)) + 
         geom_point()

ggplot(HKdatwI2021NoBP, aes(y = Average_bed_depth_cm, x = locexpdepth)) + 
  geom_boxplot()

ggplot(HKdatwI2021NoBP, aes(y = Average_bed_depth_cm, x = sitedepth)) + 
  geom_boxplot()

ggplot(HKdatwI2021NoBP, aes(y = Widest_span_cm2, x = locexpdepth)) + 
  geom_boxplot()

lengthbeddepth <- glmmTMB(log(Average_bed_depth_cm) ~ Widest_span_cm2 + 
                           (1|locexpdepth/sitedepth),
                          family = gaussian, data = HKdatwI2021NoBP)

lengthbeddepth2 <- glmmTMB(log(Average_bed_depth_cm) ~ Widest_span_cm2 + 
                             locexpdepth + 
                            locexpdepth:Widest_span_cm2 + (1|sitedepth),
                    family = gaussian, data = HKdatwI2021NoBP)

car::Anova(lengthbeddepth) # If you add in interaction, locexpdepth:Widest_span_cm2, it is sig. 
# The relationship between rubble size and bed thickness doesn't hold true in places like the reef flat.
# And it also wouldn't in blue pools.

r.squaredGLMM(lengthbeddepth) # the amount explained by the fixed effects only is very low.

# When Widest span is 0cm, the bed thickness is exp(2.375762) (10.8cm)
# This is not significantly different from 0 (because P value for intercept is <0.05)
# Then, for every 1 unit increase in Widest span, the max gap span increases by exp(0.015243) (1 cm)

# Plot

lengthbeddepth.pred <- emmeans(lengthbeddepth, ~ Widest_span_cm2,
                               at=lengthbranch.grid, type = 'response') %>% as.data.frame 

(lengthbeddepth.plot <- lengthbeddepth.pred %>% ggplot(aes(x = Widest_span_cm2, y = response)) +
    geom_point(aes(x = Widest_span_cm2, y = Average_bed_depth_cm), data = HKdatwI2021NoBP, colour = "black",
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=lower.CL, ymax =upper.CL), alpha =0.3) +
    geom_line(aes()) +
    theme_classic() + 
    xlab(bquote('Rubble length (cm)')) +
    ylab(bquote('Rubble bed thickness (cm)')) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))

# Explore differences across habitats

HKdatwI2021NoBP %>% ggplot(aes(x = Widest_span_cm2, y = Average_bed_depth_cm)) +
  geom_point(aes(x = Widest_span_cm2, y = Average_bed_depth_cm, colour = locexpdepth)) + # raw data points
  geom_smooth(aes(x = Widest_span_cm2, y = Average_bed_depth_cm, colour = locexpdepth),
              method = "lm") + # raw data points
  theme_classic() + 
  xlab(bquote('Rubble length (cm)')) +
  ylab(bquote('Rubble bed thickness (cm)')) +
  theme_minimal_hgrid(line_size = 0.2) +
  theme(text = element_text(size=10, face="bold")) + 
  theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0))

# Will remodel as reef slope vs reef flat for ease.

car::Anova(lengthbeddepth2) 

length.grid3 <- with(HKdatwI2021NoBP, 
                     list(Widest_span_cm2 = c(5, 10, 20, 30)),
                     locexpdepth = levels(locexpdepth))

emmeans(lengthbeddepth2, pairwise ~ Widest_span_cm2 | locexpdepth, at = length.grid3, 
        type = "response")

summary(lengthbeddepth2)

HKdatwI2021NoBP <- HKdatwI2021NoBP %>%
  mutate(location_type = ifelse(grepl("flat", locexpdepth), "Reef flat", "Reef slope"))

HKdatwI2021NoBP$location_type <- as.factor(HKdatwI2021NoBP$location_type)

# Plot this:

HKdatwI2021NoBP %>% ggplot(aes(x = Widest_span_cm2, y = Average_bed_depth_cm)) +
  geom_point(aes(x = Widest_span_cm2, y = Average_bed_depth_cm, colour = location_type)) + # raw data points
  geom_smooth(aes(x = Widest_span_cm2, y = Average_bed_depth_cm, colour = location_type),
              method = "lm") + # raw data points
  theme_classic() + 
  xlab(bquote('Rubble length (cm)')) +
  ylab(bquote('Rubble bed thickness (cm)')) +
  theme_minimal_hgrid(line_size = 0.2) +
  theme(text = element_text(size=10, face="bold")) + 
  theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0))


lengthbeddepth3 <- glmmTMB(log(Average_bed_depth_cm) ~ Widest_span_cm2 * location_type + (1|locexpdepth/sitedepth),
                           family = gaussian, data = HKdatwI2021NoBP)

car::Anova(lengthbeddepth3) # sig. interaction

length.grid4 <- with(HKdatwI2021NoBP, 
                     list(Widest_span_cm2 = c(5, 10, 15, 20, 25, 30)),
                     location_type = levels(location_type))

lengthbeddepth3comps <- emmeans(lengthbeddepth3, pairwise ~ Widest_span_cm2 | location_type, at = length.grid4, 
        type = "response")

save_tables_to_docx(lengthbeddepth3comps, "lengthbeddepth3comps_emmeans.docx", "lengthbeddepth3comps_comps.docx")

# The reef flat bed thickness does not increase as rubble size increases.
# Whereas it does on the reef slope.

summary(lengthbeddepth3)

r.squaredGLMM(lengthbeddepth3)

# Plot

lengthbeddepth3.grid <- with(HKdatwI2021NoBP, 
                    list(Widest_span_cm2 = seq(min(Widest_span_cm2, na.rm = TRUE), 
                                               max(Widest_span_cm2, na.rm = TRUE),len = 100),
                    location_type = levels(location_type)))


lengthbeddepth3.pred1 <- emmeans(lengthbeddepth3, ~  location_type  | Widest_span_cm2,
                               at=lengthbeddepth3.grid, type = 'response')

lengthbeddepth3.pred2 <-as.data.frame(lengthbeddepth3.pred)
View(lengthbeddepth3.pred2)

HKdatwI2021NoBP %>% group_by(location_type) %>%
  summarise(max = max(Widest_span_cm2, na.rm = TRUE), 
            min = min(Widest_span_cm2, na.rm = TRUE))

chunk1 <- lengthbeddepth3.pred2 %>%
  dplyr::group_by(location_type) %>% 
  filter((location_type=="Reef flat" & Widest_span_cm2 < 34))

chunk2 = lengthbeddepth3.pred2 %>%
  group_by(location_type) %>%
  filter((location_type=="Reef slope" & Widest_span_cm2 < 58))

colsTh <- c("Reef slope" = "black",
           "Reef flat" = "#AEE2F5")
          
pd <- position_dodge(0.5)
legend_titleTh <- "Reef zone"

(lengthbeddepth3.plot <- lengthbeddepth3.pred2 %>% ggplot(aes(x = Widest_span_cm2, y = response)) +
    geom_point(aes(x = Widest_span_cm2, y = Average_bed_depth_cm, colour = location_type), data = HKdatwI2021NoBP,
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.3) + # raw data points
    geom_ribbon(aes(ymin=lower.CL, ymax =upper.CL, fill = location_type), alpha =0.3, 
                data = chunk1, show.legend = FALSE) +
    geom_line(aes(x = Widest_span_cm2, y = response, colour = location_type), 
              data = chunk1, show.legend = FALSE) +
    geom_ribbon(aes(ymin=lower.CL, ymax =upper.CL, fill = location_type), alpha =0.3, data = chunk2) +
    geom_line(aes(x = Widest_span_cm2, y = response, colour = location_type), data = chunk2) +
    theme_classic() +
   scale_colour_manual(legend_titleTh, values=colsTh, 
                       labels=c("Reef flat", "Reef slope")) +
    scale_fill_manual(legend_titleTh, values=colsTh, 
                        labels=c("Reef flat", "Reef slope")) +
    xlab(bquote('Rubble length (cm)')) +
    ylab(bquote('Rubble bed thickness (cm)')) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="bottom") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# 4.5 x 3.3 inch

# 5.2 Branchiness ----

lengthbranch2 <- glmmTMB(Num_branches2 ~ 
                           Widest_span_cm2 * locexpdepth + (1|sitedepth/Quadrat2), 
                         family = nbinom1(), data = HKdatwI2021NoBP) # no interaction

lengthbranch <- glmmTMB(Num_branches2 ~ 
                    Widest_span_cm2 + (1|locexpdepth/sitedepth/Quadrat2), 
                  family = nbinom1(), data = HKdatwI2021NoBP)

car::Anova(lengthbranch2)

r.squaredGLMM(lengthbranch) # pretty low, 17% fixed effects, 39% fixed and random

summary(lengthbranch)

lengthbranch.grid <- with(HKdatwI2021NoBP, 
                        list(Widest_span_cm2 = seq(min(Widest_span_cm2, na.rm = TRUE), 
                                                   max(Widest_span_cm2, na.rm = TRUE),len = 100)))

lengthbranch.grid2 <- with(HKdatwI2021NoBP, 
                         list(Widest_span_cm2 = c(3, 5, 10, 15, 20, 25, 30, 35, 40)))

emmeans(lengthbranch, ~ Widest_span_cm2,
        at=lengthbranch.grid2, type = 'response')

lengthbranch.pred <- emmeans(lengthbranch, ~ Widest_span_cm2,
                           at=lengthbranch.grid, type = 'response')
lengthbranch.pred <- as.data.frame(lengthbranch.pred)

lengthbranch.pred2 <- emmeans(lengthbranch, pairwise ~ Widest_span_cm2,
                             at=lengthbranch.grid2, type = 'response')

save_tables_to_docx(lengthbranch.pred2, "lengthbranch.pred2_emmeans.docx", "lengthbranch.pred2_comps.docx")

# When Widest span is 0cm, the number of branches is exp(0.371883) (1.5)
# Then, for every 1 unit increase in Widest span, the max gap span increases by exp(0.037621) (1.04 cm)

# Plot

(sizebr1_BP.plot <- lengthbranch.pred %>% ggplot(aes(x = Widest_span_cm2, y = response)) +
    geom_point(aes(x = Widest_span_cm2, y = Num_branches2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=asymp.LCL, ymax =asymp.UCL), alpha =0.3) +
    geom_line(aes()) +
    theme_classic() + 
    scale_y_continuous(breaks = c(0,4,8,12,16), 
                     limits = c(0,17)) +
    #scale_x_continuous(limits = c(0,0.6)) +
    xlab(bquote('Rubble length (cm)')) +
    ylab(bquote('Number of branches')) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# 5 x 3.7

HKdatwI2021NoBP %>% ggplot(aes(x = Widest_span_cm2, y = Num_branches2)) +
  geom_point(aes(x = Widest_span_cm2, y = Num_branches2, colour = locexpdepth)) + # raw data points
  geom_smooth(aes(x = Widest_span_cm2, y = Num_branches2, colour = locexpdepth),
              method = "lm") + # raw data points
  theme_classic() + 
  xlab(bquote('Rubble length (cm)')) +
  ylab(bquote('Number of branches')) +
  theme_minimal_hgrid(line_size = 0.2) +
  theme(text = element_text(size=10, face="bold")) + 
  theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0))

# 5.3 Void size ----

lengthvoid2 <- glmmTMB(log(Max_gap_size_to_right_cm2) ~ 
                         Widest_span_cm2 * locexpdepth + (1|sitedepth/Quadrat2),
                       family = gaussian, data = HKdatwI2021NoBP) # no interaction, no effect of habitat

lengthvoid <- glmmTMB(log(Max_gap_size_to_right_cm2) ~ 
                     Widest_span_cm2 + (1|locexpdepth/sitedepth/Quadrat2),
                   family = gaussian, data = HKdatwI2021NoBP)

car::Anova(lengthvoid2)

r.squaredGLMM(lengthvoid) # pretty low, 20% fixed effects, 33% fixed and random

lengthvoidsumm <- summary(lengthvoid)

summary(lengthvoid)

install.packages("gtsummary")
library(gtsummary)

tbl_regression(lengthvoid, exponentiate = TRUE)

# When Widest span is 0cm, the max gap span is exp(-0.002761) (1cm) (so, means it didn't get much smaller than this)
# But this is not significantly different from 0 (because P value for intercept is not <0.05), so we can say that when Widest span is 0cm, 
# the max gap also is not significantly different from 0.
# Then, for every 1 unit increase in Widest span, the max gap span increases by exp(0.038222) (1.04 cm)

lengthvoid.pred <- emmeans(lengthvoid, ~ Widest_span_cm2,
                             at=lengthbranch.grid, type = 'response') %>% as.data.frame 

# When Widest span is 0cm, the number of branches is exp(0.371883) (1.5)
# Then, for every 1 unit increase in Widest span, the max gap span increases by exp(0.037621) (1.04 cm)

# Plot

(lengthvoid.plot <- lengthvoid.pred %>% ggplot(aes(x = Widest_span_cm2, y = response)) +
    geom_point(aes(x = Widest_span_cm2, y = Max_gap_size_to_right_cm2), data = HKdatwI2021NoBP, colour ="black", 
               position=position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data points
    geom_ribbon(aes(ymin=lower.CL, ymax =upper.CL), alpha =0.3) +
    geom_line(aes()) +
    theme_classic() + 
    #scale_x_continuous(limits = c(0,0.6)) +
    xlab(bquote('Rubble length (cm)')) +
    ylab(bquote('Void size (cm)')) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="none") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0)))
# 5 x 3.7

HKdatwI2021NoBP %>% ggplot(aes(x = Widest_span_cm2, y = Num_branches2)) +
  geom_point(aes(x = Widest_span_cm2, y = Max_gap_size_to_right_cm2, colour = locexpdepth)) + # raw data points
  geom_smooth(aes(x = Widest_span_cm2, y = Max_gap_size_to_right_cm2, colour = locexpdepth),
              method = "lm") + # raw data points
  theme_classic() + 
  xlab(bquote('Rubble length (cm)')) +
  ylab(bquote('Number of branches')) +
  theme_minimal_hgrid(line_size = 0.2) +
  theme(text = element_text(size=10, face="bold")) + 
  theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.4,angle = 0))

# **Quadrat angle** -----

# Plot quadrat angle (we didn't measure the angle but recorded as flat or sloped)

colsA <- c("flat" = "#4F4F4F",
           "sloped" = "#999999")

legend_titleA <- "Slope angle" 

(angleplot <- HKdatwI2021NoI %>% filter(Quadrat_angle2 != "NA") %>%
    ggplot(aes(x=locexpdepth, fill=Quadrat_angle2)) + 
    geom_bar(position = "fill") +
    labs(x="Habitat", y = "Proportion of quadrats") +
    scale_fill_manual(legend_titleA, values = colsA, 
                      labels=c("Flat", "Sloped")) +
    scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    theme_minimal_hgrid(line_size = 0.2) +
    theme(text = element_text(size=10, face="bold")) + 
    theme(axis.text.y.left = element_text(size=10, colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position="right") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 1,hjust = 0.8,angle = 35)))
# Saved as 5.5 * 3.7


# **Benthic cover** ----

# Load data 

datB <- read.csv("Benthic_Cover_Binding_Survey_Data.csv", header =TRUE)

# Data formatting

cols1 <- c("location", "exposure", "site", "depth", "quadrat_number", "year", "surveyed_both_years",
           "quadrat_number2.highlight.ones.that.needto.be.checked..weren.t.redone.in.2023.") 

datB <- datB %>% mutate_at(cols1, factor)

# Create an exposure depth variable

datB$exposuredepth <- paste(datB$exposure, datB$depth, sep = "_")
datB$exposuredepth <- as.factor(datB$exposuredepth)

# Create a location exposure depth variable

datB$locexpdepth <- paste(datB$location, datB$exposuredepth, sep = "_")
datB$locexpdepth <- as.factor(datB$locexpdepth)

# Format levels

datB$locexpdepth <- factor(datB$locexpdepth, levels = c("Keppels_Sheltered_Shallow",
                                                        "Keppels_Exposed_Shallow",
                                                        "Heron_Sheltered_Reef Flat",
                                                        "Heron_Sheltered_Shallow",
                                                        "Heron_Sheltered_Deep",
                                                        "Heron_Intermediate_Shallow",
                                                        "Heron_Intermediate_Deep",
                                                        "Heron_Exposed_Shallow",
                                                        "Heron_Exposed_Deep",
                                                        "Heron_Intermediate_Blue  Pools"))

datB$site <- factor(datB$site, levels = c("Halfway Sheltered", "Humpy Sheltered", "Clam Bay",
                                          "Halfway Exposed", "Humpy Exposed", "Red Beach",
                                          "Coral gardens", "Coral Gardens", "Halfway", "First Point",
                                          "Eco 1", "Eco 1 ", "Eco 2", "Blue Pools"),
                    labels = c("Halfway Sheltered", "Humpy Sheltered", "Clam Bay",
                               "Halfway Exposed", "Humpy Exposed", "Red Beach",
                               "Coral Gardens", "Coral Gardens", "Halfway", "First Point",
                               "Eco 1", "Eco 1", "Eco 2", "Blue Pools"))

datB$location <- factor(datB$location, levels = c("Keppels",
                                                  "Heron"),
                        labels = c("Keppels (inshore)", 
                                   "Heron (offshore)"))
names(datB)

# Combine categories:

datB <-  datB %>% mutate(hard_coral = rowSums(across(c(branching_coral,massive_coral, 
                                                       encrusting_coral,foliose_coral,plate_coral,solitary_coral))))

datB <-  datB %>% mutate(soft_inverts = rowSums(across(c(sponge, ascidian,soft_coral, other))))

datB <-  datB %>% mutate(sand_silt = rowSums(across(c(sand, silt))))

# Save datasheet for Primer analysis

datBPrimer <- datB %>% filter(locexpdepth != "Heron_Intermediate_Shallow") %>%
  filter(locexpdepth != "Heron_Intermediate_Deep") %>%
  filter(locexpdepth != "Blue Pools_BP") %>%
  filter(year != "2023") %>% droplevels() 

write.csv(datBPrimer, "datBPrimer.csv")

datBlong <- gather(datB, key = "benthiccat", value = "percent", 
                   hard_coral,
                   soft_inverts,
                   macroalgae, 
                   rock_w_turf_new.need_to_update_2023ones, 
                   rock_w_cca_new.need_to_update_2023ones, 
                   rubble_new.need_to_update_2023ones,
                   sand_silt,
                   standing_dead_coral, na.rm = TRUE)

datBlong$benthiccat <- as.factor(datBlong$benthiccat)

datBlong$benthiccat <- factor(datBlong$benthiccat,
                              levels = c('macroalgae', 'rock_w_turf_new.need_to_update_2023ones', 'rock_w_cca_new.need_to_update_2023ones', 'soft_inverts',
                                         'hard_coral', 'standing_dead_coral','sand_silt', 'rubble_new.need_to_update_2023ones'), 
                              labels = c('Macroalgae', 'Rock with Turf', 'Rock with CCA', 'Soft inverts',
                                         'Live hard coral', 'Dead hard coral','Sand', 'Rubble'))

datBlong21 <- datBlong %>% filter(locexpdepth != "Heron_Intermediate_Shallow") %>%
  filter(locexpdepth != "Heron_Intermediate_Deep") %>%
  filter(year != "2023") %>% droplesevels() 

datBlong2SUM <- datBlong21 %>% dplyr::group_by(year, locexpdepth, benthiccat) %>%
  summarise(prop = mean(percent),
            propsd = sd(percent),
            propse = propsd / sqrt(n())) %>% as.data.frame()

typecols <- c("Macroalgae" = "#79B07A",
              "Rock with Turf" = "#B4EEB4",
              "Rock with CCA" =  "#EEB4B4",
              "Soft inverts" = "#C799C6",
              "Live hard coral" = "peachpuff3",
              "Dead hard coral" = "#919DBD",
              "Rubble" = "#8B4500",
              "Sand" = "#F0E68C")

# **Figure 3** ----

# Need to create a dataframe WITH Blue Pools, to plot the proportion of binder groups"
# It is this one: HKdatwI2021NoI

str(HKdatwI2021NoI)
levels(HKdatwI2021NoI$sitedepth)

# Binders plot ----- 

HKdatwI2021NoIlong <- gather(HKdatwI2021NoI, key = "bindernumber", value = "type", 
                             Binder1.broad_cat, 
                             Binder2.broad_cat, 
                             Binder3.broad_cat, 
                             Binder4.broad_cat, 
                             Binder5.broad_cat, 
                             Binder6.broad_cat, 
                             Binder7.broad_cat, 
                             Binder8.broad_cat, 
                             Binder9.broad_cat,
                             Binder10.broad_cat, 
                             Binder11.broad_cat, 
                             Binder12.broad_cat, 
                             Binder13.broad_cat, 
                             Binder14.broad_cat, 
                             Binder15.broad_cat, na.rm = TRUE) # removing NAs eg where there was 1 binder and NA for 2-15.

view(HKdatwI2021NoIlong)

 HKdatwI2021NoIlongY <- HKdatwI2021NoIlong %>% filter(Bound == 1) %>% droplevels() # only the bounds

class(HKdatwI2021NoIlongY$Bound)

HKdatwI2021NoIlongY$type <- as.factor(HKdatwI2021NoIlongY$type)
levels(HKdatwI2021NoIlongY$type)

HKdatwI2021NoIlongY$type <- factor(HKdatwI2021NoIlongY$type,
                                  levels = c('Cyanobacteria', 'Turf algae', 'Lobophora', 'Other macroalgae',
                                             'CCA', 'Peyssonnelia', 'Sponge',
                                             'Ascidian', 'Tunicate', 'Bryozoan', 'Hard coral',
                                             'Cemented', 'Solid'), 
                                  labels = c('Cyanobacteria', 'Turf algae', 'Lobophora', 'Other macroalgae',
                                             'Coralline algae', 'Peyssonnelia', 'Sponge',
                                             'Ascidian', 'Ascidian', 'Bryozoan', 'Hard coral',
                                             'Solid', 'Solid')) # combining tunicate and asicidian

typecols <- c("Cyanobacteria" = "#5F9EA0",
              "Turf algae" = "#B4EEB4",
              "Lobophora" = "#B8860B",
              "Other macroalgae" = "#79B07A",
              "Coralline algae" =  "#EEB4B4",
              "Peyssonnelia" = "#F5BC49",
              "Sponge" = "#F0E68C",     
              "Ascidian" = "#C799C6",
              "Bryozoan" =   "lightblue2",
              "Hard coral" = "peachpuff3",
              "Solid" = "#919DBD")

(typebar_BP <- ggplot(HKdatwI2021NoIlongY, aes(x=locexpdepth)) +
    geom_bar(aes(fill = type), position="fill", colour="black") + 
    labs(x="Habitat", y = "Proportion binder type") + 
    scale_fill_manual(name = "Binder type", values= typecols) +
    scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    theme_classic() +
    theme(text = element_text(face="bold")) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "bottom") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 0.6,hjust = 0.5,angle = 45)))

# Plot with total count

ggplot(HKdatwI2021NoIlongY, aes(x=locexpdepth)) +
  geom_bar(aes(fill = type), position="stack") + 
  labs(x="Habitat", y = "Number of binds per binder type") + 
  scale_fill_manual(name = "Binder type", values= typecols) +
  scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                              "Off Shelt Shal", "Off Shelt Deep", 
                              "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
  theme_classic() +
  theme(text = element_text(face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "bottom") +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dashed")) +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 0.6,hjust = 0.5,angle = 45))

# Plot with # binds per type per rubble piece

library(dplyr)
library(ggplot2)

# Summarize the number of binds per Unique Piece for each type and habitat
sumdat <- HKdatwI2021NoIlong %>%
  group_by(locexpdepth, type, Unique_piece) %>%
  summarise(binds_per_piece = n())

head(sumdat)

sumdat2 <- sumdat %>%
  group_by(locexpdepth, type) %>%
  summarise(avg_binds = mean(binds_per_piece))

view(sumdat)
view(sumdat2)

ggplot(sumdat2, aes(x = locexpdepth, y = avg_binds, fill = type)) +
  geom_bar(stat="identity", position="dodge") + 
  labs(x="Habitat", y = "Average number of binds per binder type") + 
  scale_fill_manual(name = "Binder type", values= typecols) +
  scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                              "Off Shelt Shal", "Off Shelt Deep", 
                              "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
  theme_classic() +
  theme(text = element_text(face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "bottom") +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dashed")) +
  theme(axis.text.x.bottom = element_text(size=10,vjust = 0.6,hjust = 0.5,angle = 45))


# View the summarized data
print(summary_data)



# PERMANOVA - analysis in Primer

# Formatting the dataframe prior to being able to run a permanova on a distance matrix
# But want to use the dataframe WITHOUT Blue Pools, using HKdatwI2021NoBP

HKdatwI2021NoBPlong <- gather(HKdatwI2021NoBP, key = "bindernumber", value = "type", 
                             Binder1.broad_cat, 
                             Binder2.broad_cat, 
                             Binder3.broad_cat, 
                             Binder4.broad_cat, 
                             Binder5.broad_cat, 
                             Binder6.broad_cat, 
                             Binder7.broad_cat, 
                             Binder8.broad_cat, 
                             Binder9.broad_cat,
                             Binder10.broad_cat, 
                             Binder11.broad_cat, 
                             Binder12.broad_cat, 
                             Binder13.broad_cat, 
                             Binder14.broad_cat, 
                             Binder15.broad_cat, na.rm = TRUE) # removing NAs eg where there was 1 binder and NA for 2-15.

HKdatwI2021NoBPlongY <- HKdatwI2021NoBPlong %>% filter(Bound == 1) %>% droplevels() # only the bounds

class(HKdatwI2021NoBPlongY$Bound)

HKdatwI2021NoBPlongY$type <- as.factor(HKdatwI2021NoBPlongY$type)
levels(HKdatwI2021NoBPlongY$type)

HKdatwI2021NoBPlongY$type <- factor(HKdatwI2021NoBPlongY$type,
                                   levels = c('Cyanobacteria', 'Turf algae', 'Lobophora', 'Other macroalgae',
                                              'CCA', 'Peyssonnelia', 'Sponge',
                                              'Ascidian', 'Tunicate', 'Bryozoan', 'Hard coral',
                                              'Cemented', 'Solid'), 
                                   labels = c('Cyanobacteria', 'Turf algae', 'Lobophora', 'Other macroalgae',
                                              'Coralline algae', 'Peyssonnelia', 'Sponge',
                                              'Ascidian', 'Ascidian', 'Bryozoan', 'Hard coral',
                                              'Solid', 'Solid')) # combining tunicate and asicidian

# We want the total per type from HKdatwI2021NoBPlongY, before we can put it back in wide form 
# We will put it into a wide form with categories of binders as the columns instead of binder_1, binder_2, binder_3 etc columns
HKdatwI2021NoBPlongY2 <- HKdatwI2021NoBPlongY %>% dplyr::group_by(locexpdepth, sitedepth, Quadrat2, Unique_piece, Widest_span_cm2,
                                                              type) %>% # 
  dplyr::summarise(typecount = n())

view(HKdatwI2021NoBPlongY2)
levels(HKdatwI2021NoBPlongY2$type)

# Now put into wide form again
HKdatwI2021NoBPwide <- HKdatwI2021NoBPlongY2 %>% spread(key = type, value = typecount)
view(HKdatwI2021NoBPwide)

# Change all the NAs to 0s
HKdatwI2021NoBPwide[is.na(HKdatwI2021NoBPwide)] <- 0

view(HKdatwI2021NoBPwide)
# Shows there were 440 pieces bound - check this

sum1 <- HKdatwI2021NoBP %>% group_by(Bound) %>% summarise(N = n())
# Yes, 440 were bound, 322 were not bound, and there are 18 NAs.

write.csv(HKdatwI2021NoBPwide, 'HKdatwI2021NoBPwide.csv') # dataframe used in Primer

# Make community and meta matrices

new.names <- make.names(names(HKdatwI2021NoBPwide))
names(HKdatwI2021NoBPwide) <- new.names

HKdatwI2021NoBPwideCOM <- HKdatwI2021NoBPwide %>% ungroup() %>% #Makes the community matrix (the response variables)
  dplyr::select(Cyanobacteria, Turf.algae, Lobophora, Other.macroalgae, 
                Coralline.algae, Peyssonnelia, Sponge, Ascidian, Bryozoan, 
                Hard.coral, Solid)

HKdatwI2021NoBPwideTRCOM <- sqrt(HKdatwI2021NoBPwideCOM)

HKdatwI2021NoBPwideMET <- HKdatwI2021NoBPwide %>% #Makes a meta matrix (the explanatory variables)
  dplyr::select(locexpdepth, sitedepth, Quadrat2, Unique_piece, Widest_span_cm2)

HKdatwI2021NoBPwideDIST <-vegdist(HKdatwI2021NoBPwideTRCOM, method = "bray")

# Homogeneity of variance check
b=betadisper(HKdatwI2021NoBPwideDIST,group=HKdatwI2021NoBPwideMET[,"locexpdepth"]) #checking homogeneity of variance of your factors. Distance matrix is first argument following by factor of your meta matrix (the factors)
anova(b)

TukeyHSD(b)

# Permanova is robust to violations of homogeneity of variance IF the design is balanced, which this is.

# Fit the model

comm.adonis <-adonis2(HKdatwI2021NoBPwideDIST ~ locexpdepth, 
                      data = HKdatwI2021NoBPwide, method = "bray", permutations = 999) %>% as.data.frame()

comm.adonis # yes, difference

compsN <- pairwiseAdonis::pairwise.adonis(HKdatwI2021NoBPwideTRCOM[,1:11],HKdatwI2021NoBPwideMET$locexpdepth) 

compsN # Keppels inshore exposed vs sheltered are different (due to more sponge at exposed and more lobophora at exposed), 
# whereas they are not different to eachother for primer analysis. 
# AND sheltered shallow IS different to sheltered deep (driven by more macroalgae at shelt shallow), exposed shallow (less sponge shelt shallow) 
# and exposed deep (more macroalgae at shelt shall).
# But reef flat is still not different to these 3 (I think because of the CCA?) except it IS different to exposed shallow.
# Due to more macroalgae and turf on reef flat and more CCA at exposed shallow.
# But, this doesn't consider random effects. Thus, analysis was then done in Primer.

write.csv(compsN, "compsN.csv")

# Simper

(simN <- with(HKdatwI2021NoBPwideMET, simper(HKdatwI2021NoBPwideTRCOM, locexpdepth, permutations = 999)))
summary(simN) # also done in Primer.


# **Fig 3 b & c** ----

# nMDS of the natural binder community

bindcomnmds <- metaMDS(HKdatwI2021NoBPwideTRCOM, distance = "bray", k=2,trymax=200, autotransform=TRUE,expand=FALSE, plot=FALSE)

bindcomnmds$stress # The stress is 0.04. A value below 0.2 is acceptable, stress is 1-R2 value.

stressplot(bindcomnmds) # Shows that our non-metric fit R2 is 0.99. There is just not a lot of data!!

bindcomnmds.sites.scores<- as.data.frame(scores(bindcomnmds, display = 'sites'))
bindcomnmds.sites.scores<- data.frame(bindcomnmds.sites.scores, HKdatwI2021NoBPwideMET)
bindcomnmds.species.scores<- as.data.frame(scores(bindcomnmds, display = 'species'))
bindcomnmds.species.scores$Binders<- row.names(bindcomnmds.species.scores)

head(bindcomnmds.sites.scores)
head(bindcomnmds.species.scores)

bindcomnmds.species.scores$Binders <- as.factor(bindcomnmds.species.scores$Binders)
levels(bindcomnmds.species.scores$Binders)

bindcomnmds.species.scores$Binders <- factor(bindcomnmds.species.scores$Binders,
                                             levels = c("Cyanobacteria",
                                                        "Turf.algae",
                                                        "Lobophora",
                                                        "Other.macroalgae",
                                                        "Coralline.algae",
                                                        "Peyssonnelia",
                                                        "Sponge",     
                                                        "Ascidian",
                                                        "Bryozoan",
                                                        "Hard.coral",
                                                        "Solid"),    
                                             labels = c("Cyanobacteria",
                                                        "Turf.algae",
                                                        "Lobophora",
                                                        "Other.macroalgae",
                                                        "Coralline.algae",
                                                        "Peyssonnelia",
                                                        "Sponge",     
                                                        "Ascidian",
                                                        "Bryozoan",
                                                        "Hard.coral",
                                                        "Solid"))

bindcomnmds.sites.scores$locexpdepth <- as.factor(bindcomnmds.sites.scores$locexpdepth)

levels(bindcomnmds.sites.scores$locexpdepth)

bindcomnmds.sites.scores$locexpdepth <- 
  factor(bindcomnmds.sites.scores$locexpdepth, 
         levels = c("Keppels (inshore)_Sheltered_Shallow","Keppels (inshore)_Exposed_Shallow",
                    "Heron (offshore)_Sheltered_Reef flat",
                    "Heron (offshore)_Sheltered_Shallow", "Heron (offshore)_Sheltered_Deep", 
                    "Heron (offshore)_Exposed_Shallow", "Heron (offshore)_Exposed_Deep", "Heron (offshore)_Intermediate_BP"),
         labels = c("In_S_Sh","In_Ex_Sh",
                    "O_S_RF",
                    "O_S_Sh", "O_S_Dp", 
                    "O_Ex_Sh", "O_Ex_Dp", "D_A"))

colsZ <- c("In_S_Sh" = "#6CA6CD",
           "In_Ex_Sh" = "#CD8C95",
           "O_S_RF" = "#AEE2F5",
           "O_S_Sh" = "#6CA6CD",
           "O_S_Dp" = "#104E8B",
           "O_Ex_Sh" = "#CD8C95",    
           "O_Ex_Dp" = "#A52A2A",
           "D_A" = "#999999")

head(bindcomnmds.sites.scores)

centroidsA <- aggregate(cbind(NMDS1, NMDS2) ~ locexpdepth, data = bindcomnmds.sites.scores, FUN = mean)
centroidsA <- data.frame(centroidsA) # pulling out just the averages which is what I will plot
legend_title <- "Habitat"

library(ggrepel)

# Habitat points ----

(plotbindcomm <-ggplot()+
   geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
   geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
   geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, color=locexpdepth, fill = locexpdepth), size = 5)+
   scale_colour_manual(legend_title, values=colsZ) +
   #geom_point(data=bindcomnmds.sites.scores, aes(y=NMDS2, x=NMDS1, color=locexpdepth, fill = locexpdepth), size = 1, alpha = 0.5)+
   #scale_colour_manual(legend_title, values=colsZ) +
   #scale_fill_manual(legend_title, values = colsZ) +
   xlim(-0.4,0.9) +
   geom_text_repel(data=centroidsA, aes(y=NMDS2, x=NMDS1, label = locexpdepth,
                                        hjust="inward",vjust="inward", color=locexpdepth), show.legend=FALSE) +
   theme_classic() +
   ylab("NMDS2") +
   xlab("NMDS1") +
   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
   geom_vline(xintercept=0, linetype="dashed") + 
   geom_hline(yintercept=0, linetype="dashed") +
   theme(text = element_text(size=15, face="bold")) + 
   theme(axis.text = element_text(colour = "black")) + 
   theme(panel.background =element_rect(colour = "black", size=1)) + 
   theme(axis.ticks.length=unit(.2,"cm")) +
   theme(legend.position = "none") )


# Binder segments plot ----

plotbindseg <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_segment(data=bindcomnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
               colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  geom_text_repel(data=bindcomnmds.species.scores, aes(y=NMDS2, x=NMDS1, 
                                                       hjust="inward",vjust="inward",
                                                       label = Binders), colour = "gray27") + 
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") # don't need the legend
plotbindseg

grid.arrange (typebar_BP, plotbindcomm, plotbindseg, nrow = 3)

# **Figure S3** ----

(typebar_benthic <- ggplot(datBlong2SUM2_2, aes(x = locexpdepth, 
                                                y = prop, fill = benthiccat)) +
    geom_bar(stat = "identity", position = "fill") + 
    labs(x="Habitat", y = "Proportion benthic cover") + 
    scale_fill_manual(name = "Benthic category", values= typecols) +
    scale_x_discrete(labels = c("In Shelt Shal","In Exp Shal",  "Off Shelt RF",
                                "Off Shelt Shal", "Off Shelt Deep", 
                                "Off Exp Shal", "Off Exp Deep", "Dep. Area")) +
    theme_classic() +
    theme(text = element_text(face="bold")) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "bottom") +
    theme(axis.text.x.bottom = element_text(size=10,vjust = 0.6,hjust = 0.5,angle = 45)))

# Summarise mean and CIs:

datBlong21SUM <- datBlong21 %>% dplyr::group_by(year, locexpdepth, benthiccat) %>%
  summarise(meanprop = mean(percent, na.rm = TRUE),
            sdprop = sd(percent, na.rm = TRUE),
            seprop = sdprop / sqrt(n()),
            ci95 = 1.96*seprop) %>%
  mutate(lowerCI = meanprop - ci95,
         upperCI = meanprop + ci95) %>% droplevels()

# Plot binding in 2021 and 2023


str(dat2)
ggplot(aes(x = locexpdepth, y = Number_binds2, colour = Trip), data = dat2) + geom_boxplot()
