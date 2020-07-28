rm(list=ls())
library(tidyverse)
library(yingtools2)
library(phyloseq)



# start here --------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(yingtools2)
library(phyloseq)
load("16s_and_metadata.RData")
#load("~/Documents/MSKCC/Projects/Simone/simone_13aug2019.RData")

#get particular color palette
get.simone.palette <- function (tax) 
{
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax.obj)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
             "Genus", "Species")
  if (!all(ranks %in% names(tax))) {
    stop("YTError: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  tax.dict <- tax[, ranks] %>% distinct()
  tax.dict$color <- rep(shades("gray", variation = 0.25), 
                        length.out = nrow(tax.dict))
  listeria <- tax.dict$Genus == "Listeria"
  tax.dict$color[listeria] <- "black"
  proteo <- tax.dict$Phylum == "Proteobacteria"
  tax.dict$color[proteo] <- rep(shades("red", variation = 0.4), 
                                length.out = sum(proteo))
  actino <- tax.dict$Phylum == "Actinobacteria"
  tax.dict$color[actino] <- rep(shades("#A77097", variation = 0.25), 
                                length.out = sum(actino))
  bacteroidetes <- tax.dict$Phylum == "Bacteroidetes"
  tax.dict$color[bacteroidetes] <- rep(shades("#51AB9B", variation = 0.25), 
                                       length.out = sum(bacteroidetes))
  clost <- tax.dict$Order == "Clostridiales"
  tax.dict$color[clost] <- rep(shades("#9C854E", variation = 0.25), 
                               length.out = sum(clost))
  lachno <- tax.dict$Family == "Lachnospiraceae"
  tax.dict$color[lachno] <- rep(shades("#EC9B96", variation = 0.25), 
                                length.out = sum(lachno))
  rumino <- tax.dict$Family == "Ruminococcaceae"
  tax.dict$color[rumino] <- rep(shades("#9AAE73", variation = 0.25), 
                                length.out = sum(rumino))
  cid.colors.new <- c(Enterococcus = "#129246", Streptococcus = "#9FB846", 
                      Staphylococcus = "#f1eb25", Lactobacillus = "#3b51a3")
  cid <- cid.colors.new[match(tax.dict$Genus, names(cid.colors.new))]
  tax.dict$color <- ifelse(is.na(cid), tax.dict$color, cid)
  tax.palette <- structure(tax.dict$color, names = as.character(tax.dict$Species))
  tax.palette
}


pal <- get.simone.palette(t %>%
                            filter(!is.na(Genus),
                                   !is.na(Species)))

# bar plots ---------------------------------------------------------------

colnames(meta_old)
unique(t$sample)
glimpse(meta_old)
meta_old <- meta_old %>%
  mutate(exp=ifelse(grepl(Sample_ID,pattern="^SB120"),"SB120","SB147"))

table(meta_old$exp)

gg <- meta_old %>%
  filter(exp=="SB120") %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group)) %>%
  filter(!is.na(mouse.number),
         !is.na(Group)) %>%
  mutate(group_sample=paste(Group,day,sep="_")) %>%
  inner_join(t) %>%
  mutate(Phylum=ifelse(Genus=="Listeria","ZZZ",Phylum)) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  group_by(group_sample) %>% 
  arrange(Species) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label = ifelse(pctseqs >= .1, as.character(Species), ""),
         tax.label=gsub(" ","\n",tax.label),
         highlight=ifelse(grepl(Genus,pattern="Listeria"),"Listeria","")) %>%
  ggplot(aes(x=reorder(group_sample,day_num),y=pctseqs,fill=Species)) +
  geom_bar(stat="identity",position="fill") +
  #scale_color_manual(values=c(NA,"Listeria"="black")) +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ Group,scales="free") +
  xlab("") +
  ggtitle("SB120") +
#  geom_text(aes(x=group_sample,y=1-y.text,label=tax.label),size=2.5,lineheight=0.6) +
  ylab("% 16S")

gg2 <- meta_old %>%
  filter(exp=="SB147") %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group)) %>%
  filter(#!is.na(mouse.number),
         !is.na(Group)) %>%
  mutate(group_sample=paste(Group,day,sep="_")) %>%
  inner_join(t) %>%
  mutate(Phylum=ifelse(Genus=="Listeria","ZZZ",Phylum)) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  group_by(group_sample) %>% 
  arrange(Species) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label = ifelse(pctseqs >= .1, as.character(Species), ""),
         tax.label=gsub(" ","\n",tax.label),
         highlight=ifelse(grepl(Genus,pattern="Listeria"),"Listeria","")) %>%
  ggplot(aes(x=reorder(group_sample,day_num),y=pctseqs,fill=Species)) +
  geom_bar(stat="identity",position="fill") +
  #scale_color_manual(values=c(NA,"Listeria"="black")) +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ Group,scales="free") +
  xlab("") +
  ggtitle("SB147") +
  #  geom_text(aes(x=group_sample,y=1-y.text,label=tax.label),size=2.5,lineheight=0.6) +
  ylab("% 16S")



library(gridExtra)
#write pdf
pdf(file="sb120_sb147_group_day_collapse_barplots_highlight_fill.pdf",
    height=10,width=10)
grid.arrange(gg,gg2)
dev.off()




# legends -----------------------------------------------------------------

#get legends... if don't want to use generic
ggl <- meta_old %>%
  filter(exp=="SB120") %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group)) %>%
  filter(!is.na(mouse.number),
         !is.na(Group)) %>%
  mutate(group_sample=paste(Group,day,sep="_")) %>%
  inner_join(t) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  group_by(group_sample) %>% 
  arrange(Species) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label = ifelse(pctseqs >= .1, as.character(Species), ""),
         tax.label=gsub(" ","\n",tax.label)) %>%
  ggplot(aes(x=reorder(group_sample,day_num),y=pctseqs,fill=Species)) +
  geom_bar(stat="identity",position="fill") +
  theme_minimal() +
  theme(#legend.position="none",
        axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ Group,scales="free") +
  xlab("") +
  ggtitle("SB120") +
  #  geom_text(aes(x=group_sample,y=1-y.text,label=tax.label),size=2.5,lineheight=0.6) +
  ylab("% 16S")


ggl2 <- meta_old %>%
  filter(exp=="SB147") %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group)) %>%
  filter(#!is.na(mouse.number),
    !is.na(Group)) %>%
  mutate(group_sample=paste(Group,day,sep="_")) %>%
  inner_join(t) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  group_by(group_sample) %>% 
  arrange(Species) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label = ifelse(pctseqs >= .1, as.character(Species), ""),
         tax.label=gsub(" ","\n",tax.label)) %>%
  ggplot(aes(x=reorder(group_sample,day_num),y=pctseqs,fill=Species)) +
  geom_bar(stat="identity",position="fill") +
  theme_minimal() +
  theme(#legend.position="none",
        axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ Group,scales="free") +
  xlab("") +
  ggtitle("SB147") +
  #  geom_text(aes(x=group_sample,y=1-y.text,label=tax.label),size=2.5,lineheight=0.6) +
  ylab("% 16S")


ggll <- erictools::extract_legend(ggl)
ggll2 <- erictools::extract_legend(ggl2)

grid.arrange(ggll)
grid.arrange(gg,gg2)

#write metadata table
if(F){
meta_old %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group)) %>%
  write.table(file="~/Documents/MSKCC/Projects/Simone/EricPlots/sb120_147_meta_data_table.txt",row.names=F,
              quote=F,sep="\t")
}



# unifrac or anonsim distance from sample 1 ------------------------------------------
library(data.table)
library(reshape2)

for_unifrac <- meta_old %>%
  filter(!is.na(Group)) %>%
  mutate(sample=paste(Sample_ID,pool,sep=".."),
         day_num=as.numeric(day),
         Group=gsub("^na.+","na_xcc_xf8ve",Group))

physub <- phy
sample_data(physub) <- set.samp(for_unifrac)
set.seed(1)
physub.r <- rarefy_even_depth(physub)
#877 sequences per samples
get.otu.melt(physub.r) %>%
  group_by(sample) %>%
  summarize(total=sum(numseqs))

uni <- distance(physub.r,method="unifrac")
uni_df <- uni %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(sample=row.names(.)) %>%
  melt(id.vars="sample") %>%
  left_join(for_unifrac) %>%
  filter(day=="-1") %>%
  select(exp_ref=exp,sample_ref=sample,day_first=day,group_ref=Group,sample=variable,unifrac=value,mouse.number_ref=mouse.number) %>%
  mutate(mouse.number_ref=as.character(mouse.number_ref),
         sample=as.character(sample),
         sample_ref=as.character(sample_ref),
         mouse.number_ref=ifelse(is.na(mouse.number_ref),str_extract(pattern="m[0-9]+",sample_ref),mouse.number_ref)) %>%
  left_join(for_unifrac %>%
              select(sample,day,Group,mouse.number,exp) %>%
              mutate(mouse.number=as.character(mouse.number),
                     sample=as.character(sample),
                     mouse.number=ifelse(is.na(mouse.number),str_extract(pattern="m[0-9]+",sample),mouse.number))) %>%
  filter(group_ref==Group,
         exp_ref==exp)
  
glimpse(uni_df)


#write unifrac distance tables
write.table(uni_df,file="unifrac_distance_table.txt",
            row.names=F,quote=F,sep="\t")

uni_up <- read.csv(file="unifrac_distance_table_updated.csv",
                   stringsAsFactors = F)
glimpse(uni_up)


uni2 <- uni_up %>%
  #mutate(mouse.number_ref=gsub("m","",mouse.number_ref),
  #       mouse.number=gsub("m","",mouse.number)) %>%
  filter(mouse.number_ref==mouse.number)


write.table(uni2, file="unifrac_distances_matched_v2.csv",
            sep=",",quote=F,row.names=F)


key <- read.delim(file="~/Documents/MSKCC/Projects/Simone/sb_key.txt") %>%
  mutate(Sample.Tube.Label=gsub("Sample ","SB.",Sample.Tube.Label),
         Sample.Tube.Label=gsub(" \\#",".",Sample.Tube.Label),
         sample_ref=paste0(Sample.Tube.Label,"..pool905")) %>%
  filter(grepl(sample_ref,pattern="147")) %>%
  select(sample_ref,mouse.number_new=mouse.number)


uni %>%
  filter(exp=="SB147",
         day_first == -1) %>%
  select(sample_ref) %>%
  unique() %>%
  left_join(key)


library(tidyverse)
library(stringr)
library(gridExtra)

uni <- read.csv("unifrac_distances_matched.csv",
                stringsAsFactors = F)

gg <- uni2 %>%
  filter(grepl(exp,pattern="120")) %>%
  mutate(day_cat=as.character(day)) %>%
  ggplot(aes(x=reorder(day_cat,day),y=unifrac)) +
  geom_point() +
  geom_line(aes(group=mouse.number)) +
  theme_minimal() +
  facet_grid(. ~ Group) +
  ggtitle("SB120") +
  xlab("Day") +
  ylab("UniFrac Distance")

gg2 <- uni2 %>%
  filter(grepl(exp,pattern="147")) %>%
  filter(!(day==-1 & unifrac > 0)) %>%
  #left_join(key)
  mutate(day_cat=as.character(day)) %>%
  ggplot(aes(x=reorder(day_cat,day),y=unifrac)) +
  geom_point() +
  geom_line(aes(group=mouse.number)) +
  theme_minimal() +
  facet_grid(. ~ Group) +
  ggtitle("SB147") +
  xlab("Day") +
  ylab("UniFrac Distance")


#write.table(gg2,file="unifrac_table_SB147_toomany.txt",row.names=F,quote=F,sep="\t")  

#make unifrac plot
grid.arrange(gg,gg2)
