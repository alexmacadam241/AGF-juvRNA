library(ragp)

data(at_nsp)

pfam_pred <- get_hmm(sequence = at_nsp$sequence[1:2],
                     id = at_nsp$Transcript.id[1:2])

pfam_pred_go <- pfam2go(data_pfam = pfam_pred, pfam = "acc")
pfam_pred_go

pfam_dtre<-read_csv("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Data/Genomes/Durusdinium_trenchii_genome/symbiodiniumPFAM2.csv")
pfam_pred_go <- pfam2go(data_pfam = pfam_dtre, pfam = "acc")
pfam_pred_go

pfam_pred_go2<- pfam_pred_go |>
  dplyr::select(id, name, acc, Pfam_name, GO_name, GO_acc)

collapsed_data <- pfam_pred_go2 %>%
  group_by(id) %>%
  summarize(GO_acc = paste(unique(GO_acc), collapse = ";"))


write.csv(collapsed_data, "/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Data/Genomes/Durusdinium_trenchii_genome/Dtre_pfam2go.csv")
