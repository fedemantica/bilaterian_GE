library(tibble)
library(dplyr)
library(tidyr)
library(MASS)
library(parallel)
library(vroom)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)
library(patchwork)

options(dplyr.summarise.inform = FALSE)

OGannotation<- read.delim("Clusted_Motif_annotation_to_OG.tab") %>% as_tibble()# File provided
OGannotation2<-OGannotation %>% dplyr::select(Filtered_OGs,Genes) %>% unique  %>% filter(Filtered_OGs !="NOT_CONSERVED_OR_EXCLUDED")%>% group_by(Filtered_OGs) %>% summarize(Genes=paste0(Genes,collapse=";"))

sizes<-vroom("All_BilaCibserved_sizes.tab",delim = "\t",col_names=c("Gene","Gene_size","Specie")) %>%
  unique() %>% dplyr::select(Gene,Specie)
# All_BilaCibserved_sizes.tab contains all Regulatory regions tested lengths.



motifmatches<-vroom("ClusteredMotifMatches_e-04.tab", delim = "\t",col_names=c("Gene","Clustered_Motif","Count_04"))
# ClusteredMotifMatches_e-04.tab is the clustred motif count per gene with a threshold of p-value<0.00001
# Scripts to obtain this file are in XXXX.bash

# Too much memory

# motifmatches_OG<-motifmatches  %>% left_join(OGannotation, by=c("Clustered_Motif")) %>% filter(!is.na(Filtered_OGs))
# rm(motifmatches)
# gc()
# motifmatches_OG<-motifmatches_OG %>% group_by(Filtered_OGs,Gene) %>% summarize(Count_04 =sum(Count_04)) %>% ungroup 
# saveRDS(motifmatches_OG."motifmatches_OG.rds")
motifmatches_OG<-readRDS("motifmatches_OG.rds") %>% split(.$Filtered_OGs)

ordered_tissues = c("Neural", "Testis", "Ovary", "Muscle", "Kidney", "Epithelial", "DigestiveTract", "Adipose")
colors_tissue<- c("Neural"="mediumblue", "Testis"="violet", "Ovary"="darkorchid",
                  "Muscle"="springgreen3", "Kidney"="gold", "Epithelial"="dodgerblue2",
                  "DigestiveTract"="chocolate1", "Adipose"="firebrick1")
v<-c("Hs2","Mm2","Bt2","Mdo","Gga","Xtr","Dre","Cmi")
i<-c("Dme","Eba","Aae", "Bmo","Tca","Ame","Bge","Cdi")

#Select only TF which tissue specifc expression is bilaterian conserved

og2check<-OGannotation %>% filter(Conservation_status=="BILATERIAN_CONSERVED") %>% dplyr::select("Filtered_OGs","Bilaterian_core") 
og2check<-og2check %>%as_tibble() %>%
  mutate(Bilaterian_core=factor(Bilaterian_core,levels=ordered_tissues)) %>%
  mutate(Test=Bilaterian_core)%>% 
  filter(!is.na(Test))%>%
  dplyr::select(Filtered_OGs,Test)%>% 
  arrange(Test) %>%  unique()
og2check_list<-split(og2check,og2check$Test)

# Loop each tissue
results<-lapply(ordered_tissues,function(x){
  if(length(og2check_list[[x]]$Filtered_OGs) ==0){return(NA);}  #Some tissues do not have TF in their core biliaterian specific genes. 
  print(x)
  directory<-paste0("/users/mirimia/fmantica/projects/bilaterian_GE/data/motif_analysis/All_version2/gene_sets_inputs/Bilateria_ancestral/",x)  #Path to where the gene lists of bilaterian tissue specific are. 
  file_list<-list.files(directory, pattern = "\\w\\w\\w-.+-.+_geneIDs.txt") #look for specific patterns from the file names
  gene2test<-lapply(file_list, function(x){
    read.table(file.path(directory,x)) %>% .$V1
  })
  allgene2test<-unlist(gene2test) #allgene2test are gene which tissue specificity is bilaterian acestral. 
  #For each TF core biliaterian specific gene we perfomed 2 different types of test, fisher exact test to identify if the proportion of tissue specific genes regulated and a linear regration.
  #We also subtest for insects and vertebrates
  TestResults<-apply(og2check_list[[x]], 1, function(line){ 
    temp<- sizes %>% left_join(motifmatches_OG[[as.character(line["Filtered_OGs"])]], by="Gene") %>%
      mutate(Count_04=as.integer(ifelse(is.na(Count_04),0,Count_04)),
             test=ifelse(Gene %in% allgene2test, "Test","bg")) 
    ################
    #### FisherTests
    ## Bilateria
    fisherRes_bila<- temp %>%
      group_by(Specie,test) %>%
      summarize(Pres_04=sum(Count_04>0),N=n(),Filtered_OGs=as.character(line["Filtered_OGs"])) %>%
      group_by(Filtered_OGs,test) %>%
      summarize(Pres_04=sum(Pres_04),N=sum(N), NOT=N-Pres_04) %>%
      dplyr::select(Filtered_OGs,test,Pres_04,NOT) %>% filter(Filtered_OGs %in% og2check_list[[x]]$Filtered_OGs) %>%
      pivot_wider(Filtered_OGs,names_from = test,values_from = c(Pres_04,NOT)) %>% rowwise() %>%
      mutate(fish_pval=fisher.test(matrix(c(Pres_04_Test,Pres_04_bg,NOT_Test,NOT_bg),ncol=2),alternative = "g")$p.value,
             N=Pres_04_Test+NOT_Test,Tissue=x,
             UnivProp=(Pres_04_Test+Pres_04_bg)/(Pres_04_Test+NOT_Test+Pres_04_bg+NOT_bg),
             ListProp=(Pres_04_Test)/(Pres_04_Test+NOT_Test)) %>% ungroup() %>% mutate(ajd=p.adjust(fish_pval,method="fdr"))
    ## Vertebrate
    fisherRes_vert<- temp %>% filter(Specie %in% v) %>%
      group_by(Specie,test) %>%
      summarize(Pres_04=sum(Count_04>0),N=n(),Filtered_OGs=as.character(line["Filtered_OGs"])) %>%
      group_by(Filtered_OGs,test) %>%
      summarize(Pres_04=sum(Pres_04),N=sum(N), NOT=N-Pres_04) %>%
      dplyr::select(Filtered_OGs,test,Pres_04,NOT) %>% filter(Filtered_OGs %in% og2check_list[[x]]$Filtered_OGs) %>%
      pivot_wider(Filtered_OGs,names_from = test,values_from = c(Pres_04,NOT)) %>% rowwise() %>%
      mutate(fish_pval=fisher.test(matrix(c(Pres_04_Test,Pres_04_bg,NOT_Test,NOT_bg),ncol=2),alternative = "g")$p.value,
             N=Pres_04_Test+NOT_Test,Tissue=x,
             UnivProp=(Pres_04_Test+Pres_04_bg)/(Pres_04_Test+NOT_Test+Pres_04_bg+NOT_bg),
             ListProp=(Pres_04_Test)/(Pres_04_Test+NOT_Test)) %>% ungroup() %>% mutate(ajd=p.adjust(fish_pval,method="fdr"))
    ## Insect
    fisherRes_inse<- temp %>% filter(Specie %in% i) %>%
      group_by(Specie,test) %>%
      summarize(Pres_04=sum(Count_04>0),N=n(),Filtered_OGs=as.character(line["Filtered_OGs"])) %>%
      group_by(Filtered_OGs,test) %>%
      summarize(Pres_04=sum(Pres_04),N=sum(N), NOT=N-Pres_04) %>%
      dplyr::select(Filtered_OGs,test,Pres_04,NOT) %>% filter(Filtered_OGs %in% og2check_list[[x]]$Filtered_OGs) %>%
      pivot_wider(Filtered_OGs,names_from = test,values_from = c(Pres_04,NOT)) %>% rowwise() %>%
      mutate(fish_pval=fisher.test(matrix(c(Pres_04_Test,Pres_04_bg,NOT_Test,NOT_bg),ncol=2),alternative = "g")$p.value,
             N=Pres_04_Test+NOT_Test,Tissue=x,
             UnivProp=(Pres_04_Test+Pres_04_bg)/(Pres_04_Test+NOT_Test+Pres_04_bg+NOT_bg),
             ListProp=(Pres_04_Test)/(Pres_04_Test+NOT_Test)) %>% ungroup() %>% mutate(ajd=p.adjust(fish_pval,method="fdr"))
    ################
    #### GLM Nbinom Test
    ## Bilateria
    temp_fit<-temp
    if(sum(temp_fit$Count_04)>0  ){
      test04<- tryCatch({MASS::glm.nb(Count_04 ~ test, data = temp_fit) %>% summary() %>% coefficients()},
                        error=function(e) NA)
    }else{test04<-NA}
    ## Vertebrate
    temp_fit<-temp%>% filter(Specie %in% v)
    if(sum(temp_fit$Count_04)>0  ){
      test04_vert<- tryCatch({MASS::glm.nb(Count_04 ~ test, data = temp_fit) %>% summary() %>% coefficients()},
                             error=function(e) NA)
    }else{test04_vert<-NA}
    ## Insect
    temp_fit<-temp%>% filter(Specie %in% i)
    if(sum(temp_fit$Count_04)>0  ){
      test04_ins<- tryCatch({MASS::glm.nb(Count_04 ~ test, data = temp_fit) %>% summary() %>% coefficients()},
                             error=function(e) NA)
    }else{test04_ins<-NA}
    
    ################
    #### Formating Output
    
    temp1<-fisherRes_bila %>% 
      dplyr::select(Filtered_OGs,Tissue,Pres_04_Test,Pres_04_bg,NOT_bg,N,fish_pval,UnivProp,ListProp) %>% 
      left_join(fisherRes_inse[,c("Filtered_OGs","Tissue","Pres_04_Test","Pres_04_bg","NOT_bg","N","fish_pval","UnivProp","ListProp")], suffix = c("",".Insects"), by=c("Filtered_OGs","Tissue")) %>%
      left_join(fisherRes_vert[,c("Filtered_OGs","Tissue","Pres_04_Test","Pres_04_bg","NOT_bg","N","fish_pval","UnivProp","ListProp")], suffix = c("",".Vertebrate"), by=c("Filtered_OGs","Tissue")) 
      
    c(temp1,
      "Beta_test"=ifelse(all(is.na(test04)),NA,test04["testTest",'Estimate']),
      "pvalue_test"=ifelse(all(is.na(test04)),NA,test04["testTest",'Pr(>|z|)']),
      "Beta_test.insect"=ifelse(all(is.na(test04_ins)),NA,test04_ins["testTest",'Estimate']),
      "pvalue_test.insect"=ifelse(all(is.na(test04_ins)),NA,test04_ins["testTest",'Pr(>|z|)']),
      "Beta_test.vertebrates"=ifelse(all(is.na(test04_vert)),NA,test04_vert["testTest",'Estimate']),
      "pvalue_test.vertebrates"=ifelse(all(is.na(test04_vert)),NA,test04_vert["testTest",'Pr(>|z|)']))
  }) %>% do.call(rbind,.)
  TestResults %>% as_tibble() %>% 
    mutate(Filtered_OGs=as.character(Filtered_OGs),Tissue=as.character(Tissue),
           UnivProp=as.numeric(UnivProp),UnivProp.Vertebrate =as.numeric(UnivProp.Vertebrate),UnivProp.Insects=as.numeric(UnivProp.Insects),
           ListProp=as.numeric(ListProp),ListProp.Vertebrate =as.numeric(ListProp.Vertebrate),ListProp.Insects=as.numeric(ListProp.Insects),
           fish_pval=as.numeric(fish_pval),fish_pval.Vertebrate=as.numeric(fish_pval.Vertebrate),fish_pval.Insects=as.numeric(fish_pval.Insects),
           Beta_test=as.numeric(Beta_test),Beta_test.vertebrates=as.numeric(Beta_test.vertebrates),Beta_test.insect=as.numeric(Beta_test.insect),
           pvalue_test=as.numeric(pvalue_test),pvalue_test.vertebrates=as.numeric(pvalue_test.vertebrates),pvalue_test.insect=as.numeric(pvalue_test.insect))
  
}) %>% do.call(rbind,.)


#Filter results

finalRes<-results %>%
  filter((pvalue_test<0.05 & Beta_test>0)  | fish_pval < 0.05) %>%
  filter((fish_pval.Insects<0.05 & fish_pval.Vertebrate <0.05) | ((pvalue_test.insect<0.05 & Beta_test.insect >0) & (pvalue_test.vertebrates<0.05  & Beta_test.vertebrates >0))) %>%
  left_join(OGannotation2, by=c("Filtered_OGs")) %>% mutate(Tissue=factor(Tissue, levels =ordered_tissues)) %>% arrange(Tissue)

#Plot results

plt1<-finalRes %>% 
  mutate(xbil=ListProp/UnivProp,
         xvert=ListProp.Vertebrate/UnivProp.Vertebrate,
         xinse=ListProp.Insects/UnivProp.Insects) %>% 
  dplyr::select(Tissue,Genes,xbil,xvert,xinse,
                fish_pval,fish_pval.Insects,fish_pval.Vertebrate,
                Beta_test,pvalue_test,
                Beta_test.insect,pvalue_test.insect,
                Beta_test.vertebrates,pvalue_test.vertebrates) %>%
  mutate(Tissue=factor(Tissue,levels = ordered_tissues),
         Genes=substr(Genes,1,10)) %>% arrange(Tissue,desc(xbil))%>%
  mutate(Genes=factor(Genes,levels = rev(Genes))) %>%
  ggplot(aes(x=Genes, color=Tissue, fill=Tissue))+
  geom_point(aes(y=xbil, size=Beta_test), shape=19) + 
  geom_point(aes(y=xvert, size=Beta_test.vertebrates), shape=6)+ 
  geom_point(aes(y=xinse, size=Beta_test.insect), shape=2)+ 
  scale_color_manual(values=colors_tissue) +  
  scale_fill_manual(values=colors_tissue) +  coord_flip() + 
  geom_hline(yintercept = 1, color="red", linetype="dashed")+
  scale_size_binned(range = c(0.5,5),breaks = c(0,0.3))+
  scale_alpha_binned(range = c(0.2,1),breaks = c(1,2,3))+
  theme_bw()+
  ylab("Found/Expected")
ggsave(plt1, filename = "Bilaterian_TF_OG.pdf",device = "pdf")




