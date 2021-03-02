# Make this an script
# Inputs = annotation table, name_of_the_column, relation/ maybe a two-column table with the old and new labels? 

# Trash only to transform my annotation table.

library("stringr")
library("tidyverse")


current_annotations <- read.table("/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/DSP_protein_annotation_Raul_20201211_odd_deleted.tsv", sep='\t', header= TRUE, row.names = 1)

str(current_annotations)

current_annotations$label <- str_replace_all(current_annotations$label,"1","normal")
current_annotations$label <- str_replace_all(current_annotations$label,"2","necrosis")
current_annotations$label <- str_replace_all(current_annotations$label,"3","epcr")
current_annotations$label <- str_replace_all(current_annotations$label,"5","crescent")
current_annotations$label <- str_replace_all(current_annotations$label,"6","fibrocellular")
current_annotations$label <- str_replace_all(current_annotations$label,"7","fibrous")

write.table( current_annotations, file="/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/DSP_protein_annotation_Raul_20201211_odd_deleted_characters_instead_numbers.tsv" , sep= "\t", quote=FALSE)

# Why the is no 4 label? 
  
current_annotations$label
str_replace_all(as.character(1:10),"1","normal")


1 = normal
1b = "odd"
2 = necrosis
3 = epcr
3b = epcr, odd
5 = crescent
6 = fibrocellular
7 = fibrous
7b = fibrous, odd

##### RNA 

current_annotations_RNA <- read.table("/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/Final_annot_RNA_20201212_odd_deleted.tsv", sep='\t', header= TRUE, row.names = 1)

current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"1","normal")
current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"2","necrosis")
current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"3","epcr")
current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"5","crescent")
current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"6","fibrocellular")
current_annotations_RNA$label <- str_replace_all(current_annotations_RNA$label,"7","fibrous")

rownames(current_annotations_RNA) <- rownames(current_annotations_RNA) %>% str_replace_all( "b_","_"  ) %>% str_replace_all( "a_","_"  )

write.table( current_annotations_RNA , file="/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/Final_annot_RNA_20201212_odd_deleted_characters_instead_numbers.tsv" , sep= "\t", quote=FALSE)
