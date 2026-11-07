## New dataset for cure fraction
## prostate cancer

setwd("/Users/grupoderisco/Desktop/PhD/papers/Cap4_MO_Gompertz/application")

dados = readRDS("datos_completos.rds") ## dataset of FOUSP

head(dados)
unique(dados$DESCTOPO)
unique(dados$M)
#unique(cancer$CIDADE_INS)
#cancer de prostata: "prostate  SOE"

# Filtering the type of cancer
# COLO DO UTERO, PROSTATE SOE, TESTICULO  SOE, COLON SIGMOIDE
# COLO DO UTERO
# COLON SOE
# PROSTATA
# 
# "COLON SOE"                 "COLON ASCENDENTE" [descartar]         "COLON SIGMOIDE"           
# "COLON DESCENDENTE" [descartar]        "COLON TRANSVERSO"[descartar]          "COLON ANGULO ESPLENICO DO"
# "COLON ANGULO HEPATICO DO"[descartar]  "COLON LESAO SOBREPOSTA DO"[descartar]

cancer = dados[(dados$DESCTOPO=="TESTICULO  SOE"),]

dim(cancer)
#head(prostate)

# Filtering covariates of interes
# covariates = c("CIDADE_INS","DTDIAG","DTTRAT","DTULTINFO","ULTINFO", "DESCTOPO",
#                "IDADE","SEXO", "CATEATEND","TRATAMENTO","RECDIST","EC",
#                "CIRURGIA","RADIO","QUIMIO","HORMONIO","TMO","IMUNO","OUTROS",
#                "CIRURANT","RADIOANT","QUIMIOANT","HORMOANT","TMOANT","IMUNOANT",
#                "CIRURAPOS","RADIOAPOS","QUIMIOAPOS","HORMOAPOS","TMOAPOS","IMUNOAPOS")
# 

cancer = cancer[,c("CIDADE_INS","DTDIAG","DTTRAT","DTULTINFO","ULTINFO", 
                   "IDADE","SEXO","EC","CIRURGIA","CIRURANT","CIRURAPOS","NAOTRAT")]

table(cancer$EC)
table(cancer$CIRURGIA)
table(cancer$CIRURANT)
table(cancer$CIRURAPOS)

cancer$time=(cancer$DTULTINFO - cancer$DTDIAG)/365.5

## Investigating the period of 2021
cancer=cancer[(cancer$DTDIAG > "2014-10-01") & (cancer$DTDIAG < "2017-05-31"),]
#cancer=cancer[(cancer$DTDIAG > "2013-01-01") & (cancer$DTDIAG < "2015-07-31"),]
#cancer=cancer[cancer$CIDADE_INS=="SAO PAULO",]
#head(cancer)
dim(cancer)

cancer$delta = 1*((cancer$ULTINFO == 3) | (cancer$ULTINFO == 4)) ## information of death
table(cancer$delta)
prop.table(table(cancer$delta))




cancer=cancer[which(!cancer$EC %in% c("X", "Y","0")), ]
unique(cancer$EC)
dim(cancer)

cancer$s1 = (cancer$EC %in% c("I", "IB","IS"))*1
cancer$s2 = (cancer$EC %in% c("IIB", "IIC", "IIA", "II"))*1
cancer$s3 = (cancer$EC %in% c("IIIC", "IIIB", "IIIA", "III"))*1


## categorizing surgery (Stage I is the reference)
cancer$surgery = ifelse(cancer$CIRURGIA==1,1,
                     ifelse(cancer$CIRURANT==1,1,ifelse(cancer$CIRURAPOS==1,1,0)))

## categorizing treatment (8 is the code for treatment done)
cancer$trat  = ifelse(cancer$NAOTRAT ==8,1,0)


## categorizing the age and clinical stage
cancer$age = ifelse(cancer$IDADE>=35,1,0)
length(cancer$age)


test_data = as.data.frame(cancer[,c("time", "delta", "s2", "s3", "trat", "age")])
test_data$time=as.numeric(test_data$time)
test_data=test_data[which(test_data[,"time"]!=0),]

dim(test_data)
head(test_data)



write.csv2(test_data, file="test_data.csv")
















