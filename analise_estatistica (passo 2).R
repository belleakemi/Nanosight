#Importar bibliotecas e scripts necessarios
library(dplyr)
library(tidyr)
library(broom)
library(rstatix)
library(DescTools)
library(jmv)
library(png)
library(grid)
library(gridExtra)
library(car)

#Definir pasta de trabalho 
setwd(seu_wd)

## Tabelas
custom_theme = ttheme_default(
  core = list(
    bg_params = list(fill = "#B0C4DE", col = "black"),  # Blue-grey background for core cells
    fg_params = list(col = "black")  # Black text for core cells
  ),
  colhead = list(
    bg_params = list(fill = "#6FA3D7", col = "black"),  # Light blue background for headers
    fg_params = list(col = "black", fontface = "bold")  # Black bold text for headers
  )
)

salvar_tabela <- function(df, png_name) {
  
  # Create the tableGrob
  table_grob <- tableGrob(df, theme = custom_theme, rows = NULL)
  
  # Calculate width and height
  width <- convertWidth(sum(table_grob$widths), "in", valueOnly = TRUE)
  height <- convertHeight(sum(table_grob$heights), "in", valueOnly = TRUE)
  
  # Save the table as a PNG image
  png(png_name, width = width, height = height, units = "in", res = 300)
  grid.draw(table_grob)
  dev.off()
}

#Importar tabelas base

nanosight_plus_sampleinfo_sem_outliers_mean <- read.csv("nanosight_plus_sampleinfo_sem_outliers_mean.csv")
nanosight_plus_sampleinfo_sem_outliers_porcentagem <- read.csv("nanosight_plus_sampleinfo_sem_outliers_porcentagem.csv")
nanosight_w1_sem_outliers_mean <- subset(nanosight_plus_sampleinfo_sem_outliers_mean, wave == "w1")
nanosight_w2_sem_outliers_mean <- subset(nanosight_plus_sampleinfo_sem_outliers_mean, wave == "w2")
nanosight_w1_sem_outliers_porcentagem <- subset(nanosight_plus_sampleinfo_sem_outliers_porcentagem, wave == "w1")
nanosight_w2_sem_outliers_porcentagem <- subset(nanosight_plus_sampleinfo_sem_outliers_porcentagem, wave == "w2")
nanosight_plus_sampleinfo_sem_outliers_concentracao <- read.csv("nanosight_plus_sampleinfo_sem_outliers_concentracao.csv")
nanosight_w1_sem_outliers_concentracao <- subset(nanosight_plus_sampleinfo_sem_outliers_concentracao, wave == "w1")
nanosight_w2_sem_outliers_concentracao <- subset(nanosight_plus_sampleinfo_sem_outliers_concentracao, wave == "w2")


#Análise de Variância

dados <- c("tamanho w1", "porcentagem EVs < 128.5 w1", "concentracao w1", "tamanho w2", "porcentagem EVs < 128.5 w2", "concentracao w2")
dados_w1 <- c("tamanho w1", "porcentagem EVs < 128.5 w1", "concentracao w1")
dados_w2 <- c("tamanho w2", "porcentagem EVs < 128.5 w2", "concentracao w2")


##Premissas da ANOVA

###Normalidade 

###Shapiro-wilk

shapiro_w1_tamanho <- shapiro_test(nanosight_w1_sem_outliers_mean$tamanho_mean_average)
shapiro_w2_tamanho <- shapiro_test(nanosight_w2_sem_outliers_mean$tamanho_mean_average)
shapiro_w1_porc <- shapiro_test(nanosight_w1_sem_outliers_porcentagem$EV_pequenas_porcentagem)
shapiro_w2_porc <- shapiro_test(nanosight_w2_sem_outliers_porcentagem$EV_pequenas_porcentagem)
shapiro_w1_conc <- shapiro_test(nanosight_w1_sem_outliers_concentracao$concentracao_real) #5.393e-09 
shapiro_w2_conc <- shapiro_test(nanosight_w2_sem_outliers_concentracao$concentracao_real) #9.005e-08
shapiro_tamanho <-  shapiro_test(nanosight_plus_sampleinfo_sem_outliers_mean$tamanho_mean_average) 
shapiro_porc <-  shapiro_test(nanosight_plus_sampleinfo_sem_outliers_porcentagem$EV_pequenas_porcentagem)
shapiro_conc <- shapiro_test(nanosight_plus_sampleinfo_sem_outliers_concentracao$concentracao_real)


shapiro_table <- rbind(shapiro_w1_tamanho, shapiro_w1_porc, shapiro_w1_conc, shapiro_w2_tamanho, shapiro_w2_porc, shapiro_w2_conc, shapiro_tamanho, shapiro_porc, shapiro_conc)
shapiro_dados <- c(dados, "shapiro_w1_w2_tamanho", "shapiro_w1_w2_porc", "shapiro_w1_w2_conc")
shapiro_table <- data.frame(shapiro_dados, shapiro_table) [-2]
shapiro_significante <- subset(shapiro_table, shapiro_table$p.value < 0.05)
salvar_tabela(shapiro_table, "shapiro.png")



####QQPlot para visualizar a distribuição das amostras em relação a uma distribuição normal

#####w1 tamanho

pdf("qqplot_mean_w1_sem_outliers.pdf")
qqnorm(nanosight_w1_sem_outliers_mean$tamanho_mean_average, main = "tamanho w1")
qqline(nanosight_w1_sem_outliers_mean$tamanho_mean_average, col = "red")
dev.off()

#####w2 tamanho 

pdf("qqplot_mean_w2_sem_outliers.pdf")
qqnorm(nanosight_w2_sem_outliers_mean$tamanho_mean_average, main = "tamanho w2")
qqline(nanosight_w2_sem_outliers_mean$tamanho_mean_average, col = "red")
dev.off()

##### w1 % vesiculas

pdf("qqplot_EV_pequenas_porcentagem_w1_sem_outliers.pdf")
qqnorm(nanosight_w1_sem_outliers_porcentagem$EV_pequenas_porcentagem, main = "% EVs < 128.5 w1")
qqline(nanosight_w1_sem_outliers_porcentagem$EV_pequenas_porcentagem, col = "red")
dev.off()

##### w2 % vesiculas

pdf("qqplot_EV_pequenas_porcentagem_w2_sem_outliers.pdf")
qqnorm(nanosight_w2_sem_outliers_porcentagem$EV_pequenas_porcentagem, main = "% EVs < 128.5 w2")
qqline(nanosight_w2_sem_outliers_porcentagem$EV_pequenas_porcentagem, col = "red")
dev.off()

##### w1 conc

pdf("qqplot_conc_w1_sem_outliers.pdf")
qqnorm(nanosight_w1_sem_outliers_concentracao$concentracao_real, main = "Concentracao w1")
qqline(nanosight_w1_sem_outliers_concentracao$concentracao_real, col = "red")
dev.off()

##### w2 conc

pdf("qqplot_conc_w2_sem_outliers.pdf")
qqnorm(nanosight_w2_sem_outliers_concentracao$concentracao_real, main = "Concentracao w2")
qqline(nanosight_w2_sem_outliers_concentracao$concentracao_real, col = "red")
dev.off()

###Homogeneidade (Teste de Levene)

levene_w1_tamanho <- leveneTest(tamanho_mean_average ~ Trajetoria, data = nanosight_w1_sem_outliers_mean)
levene_w2_tamanho <- leveneTest(tamanho_mean_average ~ Trajetoria, data = nanosight_w2_sem_outliers_mean)
levene_w1_porc <- leveneTest(EV_pequenas_porcentagem ~ Trajetoria, data = nanosight_w1_sem_outliers_porcentagem)
levene_w2_porc <- leveneTest(EV_pequenas_porcentagem ~ Trajetoria, data = nanosight_w2_sem_outliers_porcentagem)
levene_conc_batch <- leveneTest(concentracao_real ~ batch_nan, data = nanosight_plus_sampleinfo_sem_outliers_concentracao) #0.007
levene_tamanho_batch <- leveneTest(tamanho_mean_average ~ batch_nan, data = nanosight_plus_sampleinfo_sem_outliers_mean) #0.002
levene_porc_batch <- leveneTest(EV_pequenas_porcentagem ~ batch_nan, data = nanosight_plus_sampleinfo_sem_outliers_porcentagem) #9.891e-05
levene_w1_conc <- leveneTest(concentracao_real ~ Trajetoria*batch_nan, data = nanosight_w1_sem_outliers_concentracao) #0.8
levene_w2_conc <- leveneTest(concentracao_real ~ Trajetoria*batch_nan, data = nanosight_w2_sem_outliers_concentracao) #0.8


levene_table <- rbind(levene_w1_tamanho, levene_w1_porc, levene_w1_conc, levene_w2_tamanho, levene_w2_porc, levene_w2_conc, levene_conc_batch, levene_tamanho_batch, levene_porc_batch)
dados_levene <- c(dados, "levene_conc_batch", "levene_tamanho_batch", "levene_porc_batch")
dados_levene <- rep(dados_levene, each = 2)
levene_table <- data.frame(dados_levene, levene_table)
rownames(levene_table) <- NULL
levene_significante <- subset(levene_table, levene_table$Pr..F. < 0.05)
salvar_tabela(levene_significante, "levene_significante.png")



##ANOVA Welch 
###w1 sem outliers

welch_w1_tamanho <- welch_anova_test(nanosight_w1_sem_outliers_mean, zscore_mean ~ Trajetoria)
welch_w1_porc <- welch_anova_test(nanosight_w1_sem_outliers_porcentagem, zscore_porcentagem ~ Trajetoria)
welch_conc_batch <-  welch_anova_test(nanosight_plus_sampleinfo_sem_outliers_concentracao, zscore_concentracao ~ batch_nan)
welch_tamanho_batch <- welch_anova_test(nanosight_plus_sampleinfo_sem_outliers_mean, zscore_mean ~ batch_nan)
welch_porc_batch <- welch_anova_test(nanosight_plus_sampleinfo_sem_outliers_porcentagem, zscore_porcentagem ~ batch_nan)

welch_table <- rbind(welch_w1_tamanho, welch_w1_porc, welch_conc_batch, welch_tamanho_batch, welch_porc_batch)
dados_welch <- c(dados_w1[-3], "welch_conc_batch", "welch_tamanho_batch", "welch_porc_batch")
welch_table <- data.frame(dados_welch, welch_table) [-2]

welch_significante <- subset(welch_table, welch_table$p < 0.05)
salvar_tabela(welch_significante, "welch_significante.png")


###Post-hoc Games-Howell

gh_w1_tamanho <- games_howell_test(nanosight_w1_sem_outliers_mean, zscore_mean ~ Trajetoria)
gh_w1_tamanho[gh_w1_tamanho == "zscore_mean"] <- "tamanho w1"
gh_w1_porc <- games_howell_test(nanosight_w1_sem_outliers_porcentagem, zscore_porcentagem ~ Trajetoria)
gh_w1_porc[gh_w1_porc == "zscore_porcentagem"] <- "porcentagem EVs < 128.5 w1"
gh_conc_batch <- games_howell_test(nanosight_plus_sampleinfo_sem_outliers_concentracao, zscore_concentracao ~ batch_nan)
gh_conc_batch[gh_conc_batch == "zscore_concentracao"] <- "concentracao"
gh_tamanho_batch <- games_howell_test(nanosight_plus_sampleinfo_sem_outliers_mean, zscore_mean ~ batch_nan)
gh_porc_batch <- games_howell_test(nanosight_plus_sampleinfo_sem_outliers_porcentagem, zscore_porcentagem ~ batch_nan)


gh_table <- rbind(gh_w1_tamanho, gh_w1_porc, gh_conc_batch, gh_tamanho_batch, gh_porc_batch)
gh_table_significante <- subset(gh_table, gh_table$p.adj.signif != "ns") [-8]
names(gh_table_significante)[names(gh_table_significante) == ".y."] <- "dados"
salvar_tabela(gh_table_significante, "gh_table_significante.png")

##ANOVA simples

anova_w1_tamanho <- anova_test(nanosight_w1_sem_outliers_mean, zscore_mean ~ Trajetoria*batch_nan) 
anova_w1_porc <- anova_test(nanosight_w1_sem_outliers_porcentagem, zscore_porcentagem ~ Trajetoria*batch_nan)
anova_w2_tamanho <- anova_test(nanosight_w2_sem_outliers_mean, zscore_mean ~ Trajetoria*batch_nan) 
anova_w2_porc <- anova_test(nanosight_w2_sem_outliers_porcentagem, zscore_porcentagem ~ Trajetoria*batch_nan)
anova_w1_conc <- anova_test(nanosight_w1_sem_outliers_concentracao, zscore_concentracao ~ Trajetoria*batch_nan) #significante 
anova_w2_conc <- anova_test(nanosight_w2_sem_outliers_concentracao, zscore_concentracao ~ Trajetoria*batch_nan) #significante

anova_table <- rbind(anova_w1_tamanho, anova_w1_porc, anova_w1_conc, anova_w2_tamanho, anova_w2_porc, anova_w2_conc)
dados_anova <- rep(dados, each=3)
anova_table <- data.frame(dados_anova, anova_table)

anova_table_significante <- subset(anova_table, anova_table$p..05 != "") [-7]



###Post-hoc Tukey concentracao w1 e w2 (porc e tamanho na w2 nao deu significante)

tukey_w1_tamanho <- tukey_hsd(nanosight_w1_sem_outliers_mean, zscore_mean ~ Trajetoria*batch_nan)
tukey_w1_tamanho <- cbind(dados = "tamanho w1", tukey_w1_tamanho)
tukey_w1_porc <- tukey_hsd(nanosight_w1_sem_outliers_porcentagem, zscore_porcentagem ~ Trajetoria*batch_nan)
tukey_w1_porc <- cbind(dados = "porcentagem w1", tukey_w1_porc)
tukey_w1_conc <- tukey_hsd(nanosight_w1_sem_outliers_concentracao, zscore_concentracao ~ Trajetoria*batch_nan)
tukey_w1_conc <- cbind(dados = "concentracao w1", tukey_w1_conc)
tukey_w2_conc <- tukey_hsd(nanosight_w2_sem_outliers_concentracao, zscore_concentracao ~ Trajetoria*batch_nan)
tukey_w2_conc <- cbind(dados = "concentracao w2", tukey_w2_conc)

tukey_table <- rbind(tukey_w1_tamanho, tukey_w1_porc, tukey_w1_conc, tukey_w2_conc)
tukey_table_significante <- subset(tukey_table, tukey_table$p.adj.signif != "ns") [-10]


##ANOVA RM

###Separar pares

nanosight_pares_tamanho <- nanosight_plus_sampleinfo_sem_outliers_mean[duplicated(nanosight_plus_sampleinfo_sem_outliers_mean$subjectid) | duplicated(nanosight_plus_sampleinfo_sem_outliers_mean$subjectid, fromLast = TRUE), ]
nanosight_pares_porc <- nanosight_plus_sampleinfo_sem_outliers_porcentagem[duplicated(nanosight_plus_sampleinfo_sem_outliers_porcentagem$subjectid) | duplicated(nanosight_plus_sampleinfo_sem_outliers_porcentagem$subjectid, fromLast = TRUE), ]
nanosight_pares_conc <- nanosight_plus_sampleinfo_sem_outliers_concentracao[duplicated(nanosight_plus_sampleinfo_sem_outliers_concentracao$subjectid) | duplicated(nanosight_plus_sampleinfo_sem_outliers_concentracao$subjectid, fromLast = TRUE), ]

###Colocar como fator

nanosight_pares_tamanho$Trajetoria <- factor(nanosight_pares_tamanho$Trajetoria, levels = c("A", "B", "C", "D"))
nanosight_pares_tamanho$wave <- factor(nanosight_pares_tamanho$wave, levels = c("w1", "w2"))
nanosight_pares_porc$Trajetoria <- factor(nanosight_pares_porc$Trajetoria, levels = c("A", "B", "C", "D"))
nanosight_pares_porc$wave <- factor(nanosight_pares_porc$wave, levels = c("w1", "w2"))
nanosight_pares_conc$Trajetoria <- factor(nanosight_pares_conc$Trajetoria, levels = c("A", "B", "C", "D"))
nanosight_pares_conc$wave <- factor(nanosight_pares_conc$wave, levels = c("w1", "w2"))
nanosight_pares_conc$batch_nan <- factor(nanosight_pares_conc$batch_nan, levels = c("J", "I", "B"))

###para % 

anova_rm_porc <- anova_test(nanosight_pares_porc, dv = zscore_porcentagem, wid = subjectid, within = wave, between = c(Trajetoria, batch_nan))
anova_rm_porc <- as.data.frame(get_anova_table(anova_rm_porc, correction = "auto"))
anova_rm_porc <- cbind("porcentagem EVs < 128.5", anova_rm_porc)
names(anova_rm_porc)[names(anova_rm_porc) == "\"porcentagem EVs < 128.5\""] <- "dados"
anova_rm_porc_significante <- subset(anova_rm_porc, anova_rm_porc$`p<.05` != "") [-7]
salvar_tabela(anova_rm_porc_significante, "anova_rm_porc_significante.png")

###para tamanho (somente ANOVA, sem montar tabela final de resultado pq n deu significante)

anova_rm_tamanho <- anova_test(nanosight_pares_tamanho, dv = zscore_mean, wid = subjectid, within = wave, between = c(Trajetoria, batch_nan))
anova_rm_tamanho <- get_anova_table(anova_rm_tamanho, correction = "auto")

##para conc (tabela com significante, unica significancia foi da interacao batch e trajetoria)

anova_rm_conc <- anova_test(nanosight_pares_conc, dv = zscore_concentracao, wid = subjectid, within = wave, between = c(Trajetoria,batch_nan))
anova_rm_conc <- as.data.frame(get_anova_table(anova_rm_conc, correction = "auto"))
anova_rm_conc <- cbind("Concentracao", anova_rm_conc)
names(anova_rm_conc)[names(anova_rm_conc) == "Concentracao"] <- "dados"
anova_rm_conc_significante <- subset(anova_rm_conc, anova_rm_conc$`p<.05` != "") [-7]


###post hoc anova rm % (anova deu significante. OBS: Para conc n fiz, pq se o efeito de interesse nao deu significante, nao vou olhar o outro) 

scheffe <- ScheffeTest(zscore_porcentagem ~ Trajetoria*wave+batch_nan, subset = FALSE, data = nanosight_pares_porc, na.action = getOption("echo")) #sem int: B-A significante; com interacao: nao significante
scheffe_traj_porc <- as.data.frame(scheffe[[1]]) 
scheffe_wave_porc <- as.data.frame(scheffe[[2]]) 
scheffe_wave_traj_porc <- as.data.frame(scheffe[[3]])
scheffe_table_porc <- rbind(scheffe_traj_porc, scheffe_wave_porc, scheffe_wave_traj_porc)
scheffe_table_porc <- cbind(rownames(scheffe_table_porc), scheffe_table_porc)
row.names(scheffe_table_porc) <- NULL
scheffe_table_porc_significante <- subset(scheffe_table_porc, scheffe_table_porc$pval < 0.05)
salvar_tabela(scheffe_table_porc_significante, "scheffe_table_porc_significante.png") 


##Npar

#### KW com post-hoc DSCF

kw_dscf_tamanho_w1 <- anovaNP(formula = tamanho_mean_average ~ Trajetoria*batch_nan, data = nanosight_w1_sem_outliers_mean, es = TRUE, pairs = TRUE)
kw_dscf_tamanho_w2 <- anovaNP(formula = tamanho_mean_average ~ Trajetoria*batch_nan, data = nanosight_w2_sem_outliers_mean, es = TRUE, pairs = TRUE)
kw_dscf_porc_w1 <- anovaNP(data = nanosight_w1_sem_outliers_porcentagem, formula = EV_pequenas_porcentagem~Trajetoria*batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_porc_w2 <- anovaNP(data = nanosight_w2_sem_outliers_porcentagem, formula = EV_pequenas_porcentagem~Trajetoria*batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_conc_batch <- anovaNP(data = nanosight_plus_sampleinfo_sem_outliers_concentracao, formula = concentracao_real~batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_tamanho_batch <- anovaNP(data = nanosight_plus_sampleinfo_sem_outliers_mean, formula = tamanho_mean_average~batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_porc_batch <- anovaNP(data = nanosight_plus_sampleinfo_sem_outliers_porcentagem, formula = EV_pequenas_porcentagem~batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_conc_w1 <- anovaNP(data = nanosight_w1_sem_outliers_concentracao, formula = concentracao_real~Trajetoria*batch_nan, es = TRUE, pairs = TRUE) 
kw_dscf_conc_w2 <- anovaNP(data = nanosight_w2_sem_outliers_concentracao, formula = concentracao_real~Trajetoria*batch_nan, es = TRUE, pairs = TRUE) 

#####KW

kw_jmv_tamanho_w1 <- as.data.frame(kw_dscf_tamanho_w1$table)
kw_jmv_tamanho_w2 <- as.data.frame(kw_dscf_tamanho_w2$table)
kw_jmv_porc_w1 <- as.data.frame(kw_dscf_porc_w1$table)
kw_jmv_porc_w2 <- as.data.frame(kw_dscf_porc_w2$table)
kw_jmv_conc_batch <- as.data.frame(kw_dscf_conc_batch$table)
kw_jmv_tamanho_batch <- as.data.frame(kw_dscf_tamanho_batch$table)
kw_jmv_porc_batch <- as.data.frame(kw_dscf_porc_batch$table)
kw_jmv_conc_w1 <- as.data.frame(kw_dscf_conc_w1$table) #significante
kw_jmv_conc_w2 <- as.data.frame(kw_dscf_conc_w2$table) #significante

kw_jmv_table <- rbind(kw_jmv_tamanho_w1, kw_jmv_porc_w1, kw_jmv_conc_w1, kw_jmv_tamanho_w2, kw_jmv_porc_w2, kw_jmv_conc_w2, kw_jmv_conc_batch, kw_jmv_tamanho_batch, kw_jmv_porc_batch)
dados_kw_jmv <- c(dados, "concentracaoxbatch", "tamanhoxbatch", "porcxbatch")
kw_jmv_table <- cbind(dados_kw_jmv, kw_jmv_table) [-2]
kw_jmv_significante <- subset(kw_jmv_table, kw_jmv_table$p < 0.05)
row.names(kw_jmv_significante) <- NULL
salvar_tabela(kw_jmv_significante, "kw_jmv_significante.png")

#####DSCF (so para kw significante)
dscf_tamanho_w1 <- as.data.frame(kw_dscf_tamanho_w1$comparisons[[1]])
dscf_tamanho_w1 <- cbind(dados = "tamanho w1", dscf_tamanho_w1)
dscf_porc_w1 <- as.data.frame(kw_dscf_porc_w1$comparisons[[1]])
dscf_porc_w1 <- cbind(dados = "porcentagem EVs < 128.5 w1", dscf_porc_w1)
dscf_conc_batch <- as.data.frame(kw_dscf_conc_batch$comparisons[[1]])
dscf_conc_batch <- cbind(dados = "concentracaoxbatch", dscf_conc_batch)
dscf_tamanho_batch <- as.data.frame(kw_dscf_tamanho_batch$comparisons[[1]])
dscf_tamanho_batch <- cbind(dados = "tamanhoxbatch", dscf_tamanho_batch)
dscf_porc_batch <- as.data.frame(kw_dscf_porc_batch$comparisons[[1]])
dscf_porc_batch <- cbind(dados = "porcentagemxbatch", dscf_porc_batch)
dscf_conc_w1 <- as.data.frame(kw_dscf_conc_w1$comparisons[[1]])
dscf_conc_w1 <- cbind(dados = "concentracao w1", dscf_conc_w1)
dscf_conc_w2 <- as.data.frame(kw_dscf_conc_w2$comparisons[[1]])
dscf_conc_w2 <- cbind(dados = "concentracao w2", dscf_conc_w2)


dscf_table <- rbind(dscf_tamanho_w1, dscf_porc_w1, dscf_conc_w1, dscf_conc_w2, dscf_conc_batch, dscf_tamanho_batch, dscf_porc_batch)
dscf_significante <- subset(dscf_table, dscf_table$p < 0.05)
row.names(dscf_significante) <- NULL
salvar_tabela(dscf_significante, "dscf_significante.png")
