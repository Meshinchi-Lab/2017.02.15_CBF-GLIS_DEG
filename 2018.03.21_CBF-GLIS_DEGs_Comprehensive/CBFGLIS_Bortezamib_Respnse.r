
#Jenny Smith 

#April 26, 2018 

#Purpose: Quick look at Bortezamib response in CBFA2T3-GLIS2

options(stringsAsFactors = FALSE)
library(ggplot2)
# library(survminer)
library(GGally)
library(survival)
library(gridExtra)
library(dplyr)
library(RColorBrewer)

setwd("~/RNA_seq_Analysis/2018.03.21_CBF-GLIS_DEGs_Comprehensive/")
CDE.1031 <- read.csv("~/reference_mapping-files/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_02.28.19.csv")
# CDE.1031 <- read.csv("H://reference_mapping-files/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_02.28.19.csv")

head(CDE.1031[,1:5])

sub <- CDE.1031 %>% 
  dplyr::filter(grepl("Yes|Int", CBFA2T3_GLIS2__RNASeqCalls)) %>% 
  mutate(CBFA2T3.GLIS2=ifelse(grepl("Int",CBFA2T3_GLIS2__RNASeqCalls), "Yes",CBFA2T3_GLIS2__RNASeqCalls)) %>% 
  dplyr::filter(CBFA2T3.GLIS2 != "Unknown") %>%
  dplyr::filter(!is.na(Treatment.Arm)) %>%
  dplyr::select(USI, Reg.=Patient.registration.number,CBFA2T3.GLIS2, Treatment.Arm,
         matches("MRD|OS|EFS|CR", ignore.case = FALSE))

table(sub$CBFA2T3.GLIS2)
head(sub)

table(sub$CBFA2T3.GLIS2, sub$Treatment.Arm)
#       Arm A Arm B
# Yes     9    13


m1 <- ggplot(data = sub, aes(y=Percent.MRD.at.end.of.course.1, x= Treatment.Arm , fill=Treatment.Arm)) +
  geom_violin(draw_quantiles=0.5) +
  geom_jitter(aes(color=Treatment.Arm)) +
  theme_classic() +
  scale_fill_brewer() +
  scale_color_manual(values = c("Arm A"="lightblue4", "Arm B"="darkblue")) +
  theme(text = element_text(size=14))

# m1

m2 <- ggplot(data = sub, aes(y=Percent.MRD.at.end.of.course.2, x= Treatment.Arm , fill=Treatment.Arm)) +
  geom_violin(draw_quantiles=0.5) +
  geom_jitter(aes(color=Treatment.Arm)) +
  theme_classic() +
  scale_fill_brewer() +
  scale_color_manual(values = c("Arm A"="lightblue4", "Arm B"="darkblue")) +
  theme(text = element_text(size=14))

# m2

c1 <- ggplot(data = sub, aes(x=CR.status.at.end.of.course.1, fill=Treatment.Arm)) +
  geom_histogram(stat="count", position = "dodge",color="black") + 
  labs(y="Number of Patients") +
  theme_classic() +
  scale_fill_brewer() +
  theme(text = element_text(size=14)) 
  # stat_summary(aes(y=Percent.MRD.at.end.of.course.1))
# c1


c2 <- ggplot(data = sub, aes(x=CR.status.at.end.of.course.2, fill=Treatment.Arm)) +
  geom_histogram(stat="count", position = "dodge",color="black") + 
  labs(y="Number of Patients") +
  theme_classic() +
  scale_fill_brewer() +
  theme(text = element_text(size=14))

# c2


e <- ggplot(data = sub, aes(x=EFS.event.type.ID, fill=Treatment.Arm)) +
  geom_histogram(stat="count", position = "dodge",color="black") + 
  labs(y="Number of Patients") +
  theme_classic() +
  scale_fill_brewer() +
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

# e

chisq.test(table(sub$Treatment.Arm, sub$CR.status.at.end.of.course.2)[,-3])
fisher.test(table(sub$Treatment.Arm, sub$CR.status.at.end.of.course.2)[,-3])


OS <- survfit(Surv(OS.time.days/365.25, Recoded.OS.ID) ~ Treatment.Arm, data=sub)
EFS <- survfit(Surv(EFS.time.days/365.25, Recoded.EFS.ID) ~ Treatment.Arm, data=sub)
okm <- ggsurv(OS, surv.col = c( "lightblue3","#3182BD"),ylab="Overall Survival",
              lty.est = 1, size.est = 1.25, cens.size = 4) +
  theme(text=element_text(size=14)) + 
  theme_classic()
# okm
ekm <- ggsurv(EFS, surv.col = c( "lightblue3","#3182BD"), ylab="Event-Free Survival",
              lty.est = 1, size.est = 1.25, cens.size = 4) +
  theme(text=element_text(size=14)) + 
  theme_classic()
# ekm

# grid.arrange(okm,ekm, c1,c2,m1,m2,e, ncol=2, nrow=4)
pl <- list(okm,ekm, c1,c2,m1,m2,e)
ml <- marrangeGrob(pl, ncol=2, nrow=4)
ml
ggsave("TARGET_AML_1031_CBFGLIS_Bortezamib.tiff", ml, device = "tiff", dpi=300, width = 10, height = 14)
