rm(list = ls())

setwd("Documents/lab/body_size_evol")

#Esse código vai pegar os dados de body mass e SVL para as espécies de 
#vertebrados presentes nas filogenias

#Pegar todos os dados de fontes diferentes e separar em diferentes colunas, 
#Depois, fazer média em uma coluna, mediana e máximo de todas as fontes.

library(phytools)
library(stringi)

#Função para colocar automaticamente os dados dos traits respectivos nos dados 
#(data) do grupo indicado. "data" - dataframe onde serão colocados os dados com 
#coluna "Scientific_name" e coluna do trait como "traitCol", "ref" - de onde 
#serão retirados os dados com coluna trait de nome "trait" a "data" e com
#coluna "Species"
getSppTraits <- function(data, ref, traitCol, trait) {
	for (i in 1:nrow(data)) {
		data[colnames(data) == traitCol][i, ] <- 
			ifelse(length(ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ]) == 0,
			       data[colnames(data) == traitCol][i,],
			       ifelse(is.na(ref[colnames(ref) == 
                          trait][ref$Species == data$Scientific_name[i], ]),
			              data[colnames(data) == traitCol][i, ], 
			              ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ]))
		}
  return(data)
}

dat <- read.csv("data/taxonomy_vert_20jul21.csv")

#Amphibia

dat_am <- dat[dat$Class == "Amphibia",]

#Oliveira et al. 2017 (AmphiBIO v.1) - Body_mass_g_1 e SVL_mm_1
ref_am_1 <- read.csv("data/amphibia/AmphiBIO_v1_20jul21.csv") #AmphiBIOv1

dat_am["SVL_mm_1"] <- dat_am["Body_mass_g_1"] <- NA

dat_am <- getSppTraits(dat_am, ref_am_1, "Body_mass_g_1", "Body_mass_g")
dat_am <- getSppTraits(dat_am, ref_am_1, "SVL_mm_1", "Body_size_mm")

write.csv(dat_am, "data/amphibia/BodySizeAmphibia_09set21.csv")

##################

#Reptilia

dat_sq <- dat[dat$Class == "Reptilia",]

#Feldman et al. 2016 - Body_mass_g_1 e SVL_mm_1
ref_sq_1 <- read.csv("data/reptilia/Feldman&al2016_LepidosaurBodySizes.csv")
colnames(ref_sq_1)[colnames(ref_sq_1) == "binomial"] <- "Species"

#Meiri 2019a (Endotherms) - Body_mass_g_2
ref_sq_2 <- read.csv("data/reptilia/Meiri2019_Entotherm_blz138_suppl_appendix_1.csv") 
colnames(ref_sq_2)[colnames(ref_sq_2) == "ï..species"] <- "Species"

#Myhrvold et al. 2015 (Amniote Database) - Body_mass_3 e SVL_mm_3
ref_sq_3 <- read.csv("data/reptilia/AmnioteDatabase2015_02set2021.csv") 
ref_sq_3["Species"] <- NA
for (i in 1:nrow(ref_sq_3)) {
  ref_sq_3$Species[i] <- paste0(ref_sq_3$genus[i], "_", ref_sq_3$species[i])
}
ref_sq_3$adult_body_mass_g[ref_sq_3$adult_body_mass_g == -999] <- NA
ref_sq_3$adult_svl_cm[ref_sq_3$adult_svl_cm == -999] <- NA
ref_sq_3$adult_svl_cm <- (ref_sq_3$adult_svl_cm)*10 #SVL em cm, transformar em mm

#Meiri 2019b (LizardData v.1) - SVL_mm_4
ref_sq_4 <- read.csv("data/reptilia/Meiri2019_AppendixS1_Lizard_Data_1.0.csv") 
colnames(ref_sq_4)[colnames(ref_sq_4) == "Binomial"] <- "Species"

dat_sq["SVL_mm_max"] <- dat_sq["SVL_mm_median"] <- dat_sq["SVL_mm_mean"] <- 
  dat_sq["SVL_mm_4"] <- dat_sq["SVL_mm_3"] <- dat_sq["SVL_mm_1"] <- 
  dat_sq["Body_mass_g_max"] <- dat_sq["Body_mass_g_median"] <- 
  dat_sq["Body_mass_g_mean"] <- dat_sq["Body_mass_g_3"] <- dat_sq["Body_mass_g_2"] <- 
  dat_sq["Body_mass_g_1"] <- NA

dat_sq <- getSppTraits(dat_sq, ref_sq_1, "Body_mass_g_1", "mass..g.")
dat_sq <- getSppTraits(dat_sq, ref_sq_1, "SVL_mm_1", "max.length..mm.")
dat_sq <- getSppTraits(dat_sq, ref_sq_2, "Body_mass_g_2", "adult.mass..g.")
dat_sq <- getSppTraits(dat_sq, ref_sq_3, "Body_mass_g_3", "adult_body_mass_g")
dat_sq <- getSppTraits(dat_sq, ref_sq_3, "SVL_mm_3", "adult_svl_cm")
dat_sq <- getSppTraits(dat_sq, ref_sq_4, "SVL_mm_4", "maximum.SVL")

for (i in 1:nrow(dat_sq)) {
  x <- c(dat_sq$Body_mass_g_1[i], dat_sq$Body_mass_g_2[i], 
         dat_sq$Body_mass_g_3[i])
  x <- as.numeric(x)
  if (length(is.na(x)[is.na(x) == "TRUE"]) != 3) {
    dat_sq$Body_mass_g_mean[i] <- mean(x, na.rm = T)
    dat_sq$Body_mass_g_median[i] <- median(x, na.rm = T)
    dat_sq$Body_mass_g_max[i] <- max(x, na.rm = T)
  } else {
    dat_sq$Body_mass_g_mean[i] <- NA
    dat_sq$Body_mass_g_median[i] <- NA
    dat_sq$Body_mass_g_max[i] <- NA
  }
  
  y <- c(dat_sq$SVL_mm_1[i], dat_sq$SVL_mm_3[i], dat_sq$SVL_mm_4[i])
  y <- as.numeric(y)
  if (length(is.na(y)[is.na(y) == "TRUE"]) != 3) {
    dat_sq$SVL_mm_mean[i] <- mean(y, na.rm = T)
    dat_sq$SVL_mm_median[i] <- median(y, na.rm = T)
    dat_sq$SVL_mm_max[i] <- max(y, na.rm = T)
  } else {
    dat_sq$SVL_mm_mean[i] <- NA
    dat_sq$SVL_mm_median[i] <- NA
    dat_sq$SVL_mm_max[i] <- NA
  }
}

write.csv(dat_sq, "data/reptilia/BodySizeReptilia_09set21.csv")

#só para dividir nos táxons normais

dat_sq <- read.csv("data/reptilia/BodySizeReptilia_09set21.csv")

taxonRe <- data.frame(Family = c("Anguidae", "Anniellidae", "Diploglossidae",
                                 "Helodermatidae", "Lanthanotidae", 
                                 "Shinisauridae", "Varanidae", "Xenosauridae", 
                                 #Anguimorpha
                                 "Amphisbaenidae", "Bipedidae", "Blanidae",
                                 "Cadeidae", "Gymnophthalmidae", "Lacertidae",
                                 "Rhineuridae", "Teiidae", "Trogonophiidae", 
                                 #Lacertoidea
                                 "Carphodactylidae", "Diplodactylidae", 
                                 "Eublepharidae", "Gekkonidae", 
                                 "Phyllodactylidae", "Pygopodidae", 
                                 "Sphaerodactylidae", #Gekkota
                                 "Cordylidae", "Gerrhosauridae", "Scincidae",
                                 "Xantusiidae", #Scincoidea
                                 "Agamidae", "Chamaeleonidae", "Corytophanidae",
                                 "Crotaphytidae", "Dactyloidae", 
                                 "Hoplocercidae", "Iguanidae", "Leiocephalidae",
                                 "Leiosauridae", "Liolaemidae", "Opluridae", 
                                 "Phrynosomatidae", "Polychrotidae", 
                                 "Tropiduridae", #Iguania
                                 "Acrochordidae", "Aniliidae", "Anomalepididae",
                                 "Anomochilidae", "Boidae", "Bolyeridae", 
                                 "Colubridae", "Cylindrophiidae", "Dipsadidae",
                                 "Elapidae", "Gerrhopilidae",  "Homalopsidae",
                                 "Lamprophiidae", "Leptotyphlopidae", 
                                 "Loxocemidae", "Natricidae", "Pareatidae", 
                                 "Pseudoxenodontidae", "Pythonidae", 
                                 "Tropidophiidae", "Typhlopidae", "Uropeltidae",
                                 "Viperidae", "Xenodermatidae", "Xenopeltidae",
                                 "Xenophidiidae", "Xenotyphlopidae" #Serpentes
                                 ),
                      Taxon = c(rep("Anguimorpha", 8), rep("Lacertoidea", 9),
                                rep("Gekkota", 7), rep("Scincoidea", 4),
                                rep("Iguania", 14), rep("Serpentes", 27)))

#"Dibamidae"*, "Gymnophthalmidae", "Teiidae"  ???

dat_sq["Taxon"] <- NA
for (i in 1:nrow(dat_sq)) {
  dat_sq$Taxon[i] <- ifelse(dat_sq$Family[i] %in% taxonRe$Family,
                         taxonRe$Taxon[taxonRe$Family == dat_sq$Family[i]], NA)
}

write.csv(dat_sq, "data/reptilia/BodySizeReptilia_15set21.csv")

##################

#Aves

dat_av <- dat[dat$Class == "Aves",]

#EltonTraits - Body_mass_g_1
ref_av_1 <- read.csv("data/aves/EltonTraitsAv_20jul21.csv")
colnames(ref_av_1)[colnames(ref_av_1) == "Scientific"] <- "Species"

#Ocampo et al. 2021 - Body_mass_g_2, Body_mass_g_M_2 e Body_mass_g_F_2
ref_av_2 <- read.csv("data/aves/Ocampo&al2020_DataS1_BodyMass.csv")

#Lislevand et al. 2007 - Body_mass_g_3, Body_mass_g_M_3 e Body_mass_g_F_3
ref_av_3 <- read.csv("data/aves/Lislevand2007_30may22.csv") 
colnames(ref_av_3)[colnames(ref_av_3) == "Species_name"] <- "Species"
ref_av_3$Species <- stri_replace_all_fixed(ref_av_3$Species, " ", "_")
ref_av_3$M_mass[ref_av_3$M_mass == -999] <- NA
ref_av_3$F_mass[ref_av_3$F_mass == -999] <- NA

ref_av_3["Mean_mass"] <- NA
for (i in 1:nrow(ref_av_3)) {
  x <- ref_av_3$M_mass[i]
  y <- ref_av_3$F_mass[i]
  if (!is.na(x) & !is.na(y)) {
    ref_av_3$Mean_mass[i] <- mean(c(x,y))
  } else {
    if (is.na(x)) {
      ref_av_3$Mean_mass[i] <- y
    } else {
      if (is.na(y)) {
        ref_av_3$Mean_mass[i] <- x
      } else {
        ref_av_3$Mean_mass[i] <- NA
      }
    }
  }
}

#Myhrvold et al. 2015 (Amniote Database) - Body_mass_g_4, Body_mass_g_M_4 e 
##Body_mass_g_F_4
ref_av_4 <- read.csv("data/aves/AmnioteDatabase2015_02set2021.csv") 
ref_av_4["Species"] <- NA
for (i in 1:nrow(ref_av_4)) {
  ref_av_4$Species[i] <- paste0(ref_av_4$genus[i], "_", ref_av_4$species[i])
}
ref_av_4$adult_body_mass_g[ref_av_4$adult_body_mass_g == -999] <- NA
ref_av_4$female_body_mass_g[ref_av_4$female_body_mass_g == -999] <- NA
ref_av_4$male_body_mass_g[ref_av_4$male_body_mass_g == -999] <- NA

dat_av["Body_mass_g_F_max"] <- dat_av["Body_mass_g_F_median"] <- 
  dat_av["Body_mass_g_F_mean"] <- dat_av["Body_mass_g_F_4"] <- 
  dat_av["Body_mass_g_F_3"] <- dat_av["Body_mass_g_F_2"] <- 
  dat_av["Body_mass_g_M_max"] <- dat_av["Body_mass_g_M_median"] <- 
  dat_av["Body_mass_g_M_mean"] <- dat_av["Body_mass_g_M_4"] <- 
  dat_av["Body_mass_g_M_3"] <- dat_av["Body_mass_g_M_2"] <- 
  dat_av["Body_mass_g_max"] <- dat_av["Body_mass_g_median"] <- 
  dat_av["Body_mass_g_mean"] <- dat_av["Body_mass_g_4"] <- dat_av["Body_mass_g_3"] <- 
  dat_av["Body_mass_g_2"] <- dat_av["Body_mass_g_1"] <- NA

dat_av <- getSppTraits(dat_av, ref_av_1, "Body_mass_g_1", "BodyMass.Value")
dat_av <- getSppTraits(dat_av, ref_av_2, "Body_mass_g_2", 
                    "Average.weight..gr...Species.")
dat_av <- getSppTraits(dat_av, ref_av_2, "Body_mass_g_M_2", 
                    "Average.weight..gr...Male.")
dat_av <- getSppTraits(dat_av, ref_av_2, "Body_mass_g_F_2", 
                    "Average.weight..gr...Female.")
dat_av <- getSppTraits(dat_av, ref_av_3, "Body_mass_g_3", "Mean_mass")
dat_av <- getSppTraits(dat_av, ref_av_3, "Body_mass_g_M_3", "M_mass")
dat_av <- getSppTraits(dat_av, ref_av_3, "Body_mass_g_F_3", "F_mass")
dat_av <- getSppTraits(dat_av, ref_av_4, "Body_mass_g_4", "adult_body_mass_g")
dat_av <- getSppTraits(dat_av, ref_av_4, "Body_mass_g_M_4", "male_body_mass_g")
dat_av <- getSppTraits(dat_av, ref_av_4, "Body_mass_g_F_4", "female_body_mass_g")

for (i in 1:nrow(dat_av)) {
  x <- c(dat_av$Body_mass_g_1[i], dat_av$Body_mass_g_2[i], dat_av$Body_mass_g_3[i],
       dat_av$Body_mass_g_4[i])
  x <- as.numeric(x)
  if (length(is.na(x)[is.na(x) == "TRUE"]) != 4) {
    dat_av$Body_mass_g_mean[i] <- mean(x, na.rm = T)
    dat_av$Body_mass_g_median[i] <- median(x, na.rm = T)
    dat_av$Body_mass_g_max[i] <- max(x, na.rm = T)
  } else {
    dat_av$Body_mass_g_mean[i] <- NA
    dat_av$Body_mass_g_median[i] <- NA
    dat_av$Body_mass_g_max[i] <- NA
  }
  
  y <- c(dat_av$Body_mass_g_M_2[i], dat_av$Body_mass_g_M_3[i], 
       dat_av$Body_mass_g_M_4[i])
  y <- as.numeric(y)
  if (length(is.na(y)[is.na(y) == "TRUE"]) != 3) {
    dat_av$Body_mass_g_M_mean[i] <- mean(y, na.rm = T)
    dat_av$Body_mass_g_M_median[i] <- median(y, na.rm = T)
    dat_av$Body_mass_g_M_max[i] <- max(y, na.rm = T)
  } else {
    dat_av$Body_mass_g_M_mean[i] <- NA
    dat_av$Body_mass_g_M_median[i] <- NA
    dat_av$Body_mass_g_M_max[i] <- NA
  }
  
  z <- c(dat_av$Body_mass_g_F_2[i], dat_av$Body_mass_g_F_3[i], 
       dat_av$Body_mass_g_F_4[i])
  z <- as.numeric(z)
  if (length(is.na(z)[is.na(z) == "TRUE"]) != 3) {
    dat_av$Body_mass_g_F_mean[i] <- mean(z, na.rm = T)
    dat_av$Body_mass_g_F_median[i] <- median(z, na.rm = T)
    dat_av$Body_mass_g_F_max[i] <- max(z, na.rm = T)
  } else {
    dat_av$Body_mass_g_F_mean[i] <- NA
    dat_av$Body_mass_g_F_median[i] <- NA
    dat_av$Body_mass_g_F_max[i] <- NA
  }
}

write.csv(dat_av, "data/aves/BodySizeAves_10set22.csv")

##################

#Mamíferos

dat_ma <- dat[dat$Class == "Mammalia",]

#EltonTraits - Body_mass_g_1
ref_ma_1 <- read.csv("data/mammalia/EltonTraitsMa_20jul21.csv") 
colnames(ref_ma_1)[colnames(ref_ma_1) == "Scientific"] <- "Species"

#Phylacine - Body_mass_g_2
ref_ma_2 <- read.csv("data/mammalia/Phylacine_Trait_data_20jul21.csv") 
colnames(ref_ma_2)[colnames(ref_ma_2) == "Binomial.1.2"] <- "Species"

#Phanteria - Body_mass_g_3
ref_ma_3 <- read.table("data/mammalia/PanTHERIA_1-0_WR93_Aug2008.txt", 
                      sep = "\t", header = T)
colnames(ref_ma_3)[colnames(ref_ma_3) == "MSW93_Binomial"] <- "Species"
ref_ma_3$Species <- stri_replace_all_fixed(ref_ma_3$Species, " ", "_")
ref_ma_3$X5.1_AdultBodyMass_g[ref_ma_3$X5.1_AdultBodyMass_g == -999] <- NA

#Ocampo et al. 2021 - Body_mass_g_4
ref_ma_4 <- read.csv("data/mammalia/Ocampo&al2020_DataS1_BodyMass.csv") 

dat_ma["Body_mass_g_max"] <- dat_ma["Body_mass_g_median"] <- 
  dat_ma["Body_mass_g_mean"] <- dat_ma["Body_mass_g_4"] <- dat_ma["Body_mass_g_3"] <- 
  dat_ma["Body_mass_g_2"] <- dat_ma["Body_mass_g_1"] <- NA

dat_ma <- getSppTraits(dat_ma, ref_ma_1, "Body_mass_g_1", "BodyMass.Value")
dat_ma <- getSppTraits(dat_ma, ref_ma_2, "Body_mass_g_2", "Mass.g")
dat_ma <- getSppTraits(dat_ma, ref_ma_3, "Body_mass_g_3", "X5.1_AdultBodyMass_g")
dat_ma <- getSppTraits(dat_ma, ref_ma_4, "Body_mass_g_4", 
                    "Average.weight..gr...Species.")

for (i in 1:nrow(dat_ma)) {
  x <- c(dat_ma$Body_mass_g_1[i], dat_ma$Body_mass_g_2[i], dat_ma$Body_mass_g_3[i],
       dat_ma$Body_mass_g_4[i])
  x <- as.numeric(x)
  if (length(is.na(x)[is.na(x) == "TRUE"])!=4) {
    dat_ma$Body_mass_g_mean[i] <- mean(x, na.rm = T)
    dat_ma$Body_mass_g_median[i] <- median(x, na.rm = T)
    dat_ma$Body_mass_g_max[i] <- max(x, na.rm = T)
  } else {
    dat_ma$Body_mass_g_mean[i] <- NA
    dat_ma$Body_mass_g_median[i] <- NA
    dat_ma$Body_mass_g_max[i] <- NA
  }
}

write.csv(dat_ma, "data/mammalia/BodySizeMammalia_09set21.csv")
