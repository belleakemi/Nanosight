library(rstatix)
library(car)
library(WRS2)
library(dplyr)
library(dunn.test)
library(grid)
library(gridExtra)
library(PMCMRplus)
library(gt)
library(webshot2)

#Definir pasta de trabalho
setwd("C:/Users/Belle/Documents/Belle - Nanosight")

###ESTATÍSTICAS NOVAS
#Qui-quadrado para variáveis categóricas
table_batch_w1 <- table(nanosight_plus_sampleinfo_equalized_w1$Trajetoria, nanosight_plus_sampleinfo_equalized_w1$batch_nan)
chisq.test(table_batch_w1) #p-value = 7.555e-07
table_batch_w2 <- table(nanosight_plus_sampleinfo_equalized_w2$Trajetoria, nanosight_plus_sampleinfo_equalized_w2$batch_nan)
chisq.test(table_batch_w2) #p-value = 0.0009936
#ANOVA para variáveis contínuas com distribuição normal
aov_result_tamanho_w1 <- aov(tamanho_mean_average ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w1)
summary(aov_result_tamanho_w1) #p 0.22
aov_result_tamanho_w2 <- aov(tamanho_mean_average ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w2)
summary(aov_result_tamanho_w2) #p 0.164
#Kruskal-Wallis para variáveis contínuas com distribuição fora da normalidade
kruskal.test(EV_pequenas_porcentagem ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w1) #p 0.003524
kruskal.test(EV_pequenas_porcentagem ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w2) #5.239e-05
kruskal.test(concentracao_average ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w1) #p 1.279e-05
kruskal.test(concentracao_average ~ batch_nan, data = nanosight_plus_sampleinfo_equalized_w2) #p 4.456e-05
##Modelagem estatística para ajustar o efeito batch
#Modelo misto se houver repetição ou estrutura hierárquica:
library(lme4)
lmer_result_mean_w1 <- lmer(tamanho_mean_average ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w1)
summary(lmer_result_mean_w1)
lmer_result_mean_w2 <- lmer(tamanho_mean_average ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w2)
summary(lmer_result_mean_w2)
lmer_result_porcentagem_w1 <- lmer(EV_pequenas_porcentagem ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w1)
summary(lmer_result_porcentagem_w1)
lmer_result_porcentagem_w2 <- lmer(EV_pequenas_porcentagem ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w2)
summary(lmer_result_porcentagem_w2)
lmer_result_concentracao_w1 <- lmer(concentracao_average ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w1)
summary(lmer_result_concentracao_w1)
lmer_result_concentracao_w2 <- lmer(concentracao_average ~ Trajetoria + (1|batch_nan), data = nanosight_plus_sampleinfo_equalized_w2)
summary(lmer_result_concentracao_w2)
#Correção do efeito batch
set.seed(123)
sample_i1 <- nanosight_batch_i1[sample(1:nrow(nanosight_batch_i1), size = 10), ]
equalized_data <- rbind(nanosight_batch_b, sample_i1, nanosight_batch_i2, nanosight_batch_j)
#Dados com equalized batch I1
nanosight_plus_sampleinfo_equalized <- rbind(nanosight_batch_b_i2_j, sample_i1)
nanosight_plus_sampleinfo_equalized_w1 <- subset(nanosight_plus_sampleinfo_equalized, wave == "w1")
nanosight_plus_sampleinfo_equalized_w2 <- subset(nanosight_plus_sampleinfo_equalized, wave == "w2")

### ANÁLISES ESTATÍSTICAS ###
#Testes de normalidade (premissa para o ANOVA)
shapiro_teste_tamanho_w1 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w1$tamanho_mean_average) #p = 0.05609
shapiro_teste_tamanho_w2 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w2$tamanho_mean_average) #p = 0.2391
shapiro_teste_porcentagem_w1 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w1$EV_pequenas_porcentagem) #p = 5.742e-07
shapiro_teste_porcentagem_w2 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w2$EV_pequenas_porcentagem) #p = 1.728e-07
shapiro_teste_concentracao_w1 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w1$concentracao_average) #p = 4.375e-10
shapiro_teste_concentracao_w2 <- shapiro.test(nanosight_plus_sampleinfo_equalized_w2$concentracao_average) #p = 4.446e-08

#Testes de homogeneidade de variâncias (premissa para o ANOVA) #dá um aviso de variável, mas o resultado dá certo
modelo_w1_mean <- aov(tamanho_mean_average~Trajetoria, data = nanosight_intersect_w1)
levene_w1_mean <- leveneTest(modelo_w1_mean) #p = 0.003557343 (não homogêneas)
modelo_w2_mean <- aov(tamanho_mean_average~Trajetoria, data = nanosight_intersect_w2)
levene_w2_mean <- leveneTest(modelo_w2_mean) #p = 0.1708796 (homogêneas)
modelo_w1_porcentagem <- aov(EV_pequenas_porcentagem~Trajetoria, data = nanosight_intersect_w1)
levene_w1_porcentagem <- leveneTest(modelo_w1_porcentagem) #p = 0.0001411936 (não homogêneas)
modelo_w2_porcentagem <- aov(EV_pequenas_porcentagem~Trajetoria, data = nanosight_intersect_w2)
levene_w2_porcentagem <- leveneTest(modelo_w2_porcentagem) #p = 0.01902797 (não homogêneas)
modelo_w1_concentracao <- aov(concentracao_real~Trajetoria, data = nanosight_intersect_w1)
levene_w1_concentracao <- leveneTest(modelo_w1_concentracao) #p = 0.2657289 (homogêneas)
modelo_w2_concentracao <- aov(concentracao_real~Trajetoria, data = nanosight_intersect_w2)
levene_w2_concentracao <- leveneTest(modelo_w2_concentracao) #p = 0.3430197 (homogêneas)

#ANOVA robusta
modelo_w1_mean_welch <- oneway.test(tamanho_mean_average ~ Trajetoria, data = nanosight_intersect_w1, var.equal = FALSE) #0.009377 (rejeito a hipótese nula)
resultado_robusto_t1way_w2_mean <- t1way(tamanho_mean_average~Trajetoria, data = nanosight_intersect_w2) #0.3600361 (hipótese nula possível)
modelo_w1_porcentagem_welch <- oneway.test(EV_pequenas_porcentagem ~ Trajetoria, data = nanosight_intersect_w1, var.equal = FALSE) #4.655e-06 (rejeito a hipótese nula)
modelo_w2_porcentagem_welch <- oneway.test(EV_pequenas_porcentagem ~ Trajetoria, data = nanosight_intersect_w2, var.equal = FALSE) #0.0007951 (rejeito a hipótese nula)
resultado_robusto_t1way_w1_concentracao <- t1way(concentracao_real~Trajetoria, data = nanosight_intersect_w1) #0.00034 (rejeito a hipótese nula)
resultado_robusto_t1way_w2_concentracao <- t1way(concentracao_real~Trajetoria, data = nanosight_intersect_w2) #0.00482 (rejeito a hipótese nula)

#### Análises corrigindo para o efeito batch
##ANOVA robusta de 2 vias
#tamanho/batch/trajetoria na w1
resultado_anova_duas_vias_tam_batch_traj_w1 <- aov(tamanho_mean_average ~ Trajetoria * batch_nan, data = nanosight_intersect_w1)
anova_summary_tam_batch_traj_w1 <- summary(resultado_anova_duas_vias_tam_batch_traj_w1)
anova_results_df_tam_batch_traj_w1 <- as.data.frame(anova_summary_tam_batch_traj_w1[[1]])
print(anova_results_df_tam_batch_traj_w1)
pdf("plots/anova_results_df_tam_batch_traj_w1.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Size in relation to Batch and Trajectory - w1", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_tam_batch_traj_w1)  # Desenha a tabela no arquivo PDF
dev.off()
#tamanho/batch/trajetoria na w2
resultado_anova_duas_vias_tam_batch_traj_w2 <- aov(tamanho_mean_average ~ Trajetoria * batch_nan, data = nanosight_intersect_w2)
anova_summary_tam_batch_traj_w2 <- summary(resultado_anova_duas_vias_tam_batch_traj_w2)
anova_results_df_tam_batch_traj_w2 <- as.data.frame(anova_summary_tam_batch_traj_w2[[1]])
print(anova_results_df_tam_batch_traj_w2)
pdf("plots/anova_results_df_tam_batch_traj_w2.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Size in relation to Batch and Trajectory - w2", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_tam_batch_traj_w2)  # Desenha a tabela no arquivo PDF
dev.off()
#porcentagem/batch/trajetoria na w1
resultado_anova_duas_vias_porc_batch_traj_w1 <- aov(EV_pequenas_porcentagem ~ Trajetoria * batch_nan, data = nanosight_intersect_w1)
anova_summary_porc_batch_traj_w1 <- summary(resultado_anova_duas_vias_porc_batch_traj_w1)
anova_results_df_porc_batch_traj_w1 <- as.data.frame(anova_summary_porc_batch_traj_w1[[1]])
print(anova_results_df_porc_batch_traj_w1)
pdf("plots/anova_results_df_porc_batch_traj_w1.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Percentage of Small EV's in relation to Batch and Trajectory - w1", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_porc_batch_traj_w1)  # Desenha a tabela no arquivo PDF
dev.off()
#porcentagem/batch/trajetoria na w2
resultado_anova_duas_vias_porc_batch_traj_w2 <- aov(EV_pequenas_porcentagem ~ Trajetoria * batch_nan, data = nanosight_intersect_w2)
anova_summary_porc_batch_traj_w2 <- summary(resultado_anova_duas_vias_porc_batch_traj_w2)
anova_results_df_porc_batch_traj_w2 <- as.data.frame(anova_summary_porc_batch_traj_w2[[1]])
print(anova_results_df_porc_batch_traj_w2)
pdf("plots/anova_results_df_porc_batch_traj_w2.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Percentage of Small EV's in relation to Batch and Trajectory - w2", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_porc_batch_traj_w2)  # Desenha a tabela no arquivo PDF
dev.off()
#concentração/batch/trajetoria na w1
resultado_anova_duas_vias_conc_batch_traj_w1 <- aov(concentracao_average ~ Trajetoria * batch_nan, data = nanosight_intersect_w1)
anova_summary_conc_batch_traj_w1 <- summary(resultado_anova_duas_vias_conc_batch_traj_w1)
anova_results_df_conc_batch_traj_w1 <- as.data.frame(anova_summary_conc_batch_traj_w1[[1]])
print(anova_results_df_conc_batch_traj_w1)
pdf("plots/anova_results_df_conc_batch_traj_w1.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Concentration in relation to Batch and Trajectory - w1", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_conc_batch_traj_w1)  # Desenha a tabela no arquivo PDF
dev.off()
#concentração/batch/trajetoria na w2
resultado_anova_duas_vias_conc_batch_traj_w2 <- aov(concentracao_average ~ Trajetoria * batch_nan, data = nanosight_intersect_w2)
anova_summary_conc_batch_traj_w2 <- summary(resultado_anova_duas_vias_conc_batch_traj_w2)
anova_results_df_conc_batch_traj_w2 <- as.data.frame(anova_summary_conc_batch_traj_w2[[1]])
print(anova_results_df_conc_batch_traj_w2)
pdf("plots/anova_results_df_conc_batch_traj_w2.pdf", width = 10, height = 6)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Robust 2-way ANOVA for Concentration in relation to Batch and Trajectory - w2", x = 0.5, y = 0.7, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(anova_results_df_conc_batch_traj_w2)  # Desenha a tabela no arquivo PDF
dev.off()

##Teste de Tukey para comparações múltiplas entre as trajetórias
#Tamanho para batch e Trajetória na w1
nanosight_intersect_w1$batch_nan <- as.factor(nanosight_intersect_w1$batch_nan)
nanosight_intersect_w1$Trajetoria <- as.factor(nanosight_intersect_w1$Trajetoria)
tukey_test_tam_batch_traj_w1 <- tukey_hsd(resultado_anova_duas_vias_tam_batch_traj_w1, "Trajetoria:batch_nan")
summary(tukey_test_tam_batch_traj_w1)
#Tamanho para batch e Trajetória na w2
nanosight_intersect_w2$batch_nan <- as.factor(nanosight_intersect_w2$batch_nan)
nanosight_intersect_w2$Trajetoria <- as.factor(nanosight_intersect_w2$Trajetoria)
tukey_test_tam_batch_traj_w2 <- tukey_hsd(resultado_anova_duas_vias_tam_batch_traj_w2, "Trajetoria:batch_nan")
summary(tukey_test_tam_batch_traj_w2)
#Porcentagem para batch e Trajetória na w1
nanosight_intersect_w1$batch_nan <- as.factor(nanosight_intersect_w1$batch_nan)
nanosight_intersect_w1$Trajetoria <- as.factor(nanosight_intersect_w1$Trajetoria)
tukey_test_porc_batch_traj_w1 <- tukey_hsd(resultado_anova_duas_vias_porc_batch_traj_w1, "Trajetoria:batch_nan")
summary(tukey_test_porc_batch_traj_w1)
#Porcentagem para batch e Trajetória na w2
nanosight_intersect_w2$batch_nan <- as.factor(nanosight_intersect_w2$batch_nan)
nanosight_intersect_w2$Trajetoria <- as.factor(nanosight_intersect_w2$Trajetoria)
tukey_test_porc_batch_traj_w2 <- tukey_hsd(resultado_anova_duas_vias_porc_batch_traj_w2, "Trajetoria:batch_nan")
summary(tukey_test_porc_batch_traj_w2)
#Concentração para batch e Trajetória na w1
nanosight_intersect_w1$batch_nan <- as.factor(nanosight_intersect_w1$batch_nan)
nanosight_intersect_w1$Trajetoria <- as.factor(nanosight_intersect_w1$Trajetoria)
tukey_test_conc_batch_traj_w1 <- tukey_hsd(resultado_anova_duas_vias_conc_batch_traj_w1, "Trajetoria:batch_nan")
summary(tukey_test_conc_batch_traj_w1)
#Concentração para batch e Trajetória na w2
nanosight_intersect_w2$batch_nan <- as.factor(nanosight_intersect_w2$batch_nan)
nanosight_intersect_w2$Trajetoria <- as.factor(nanosight_intersect_w2$Trajetoria)
tukey_test_conc_batch_traj_w2 <- tukey_hsd(resultado_anova_duas_vias_conc_batch_traj_w2, "Trajetoria:batch_nan")
summary(tukey_test_conc_batch_traj_w2)
###Tabela dos resultados dos testes de Tukey
# Adicionar uma coluna de origem em cada tabela antes de combinar
tukey_test_tam_batch_traj_w1$Source <- "Tamanho - w1"
tukey_test_tam_batch_traj_w2$Source <- "Tamanho - w2"
tukey_test_porc_batch_traj_w1$Source <- "Porcentagem - w1"
tukey_test_porc_batch_traj_w2$Source <- "Porcentagem - w2"
tukey_test_conc_batch_traj_w1$Source <- "Concentração - w1"
tukey_test_conc_batch_traj_w2$Source <- "Concentração - w2"
#Tamanho na w1
tukey_table_tamanho_W1 <- tukey_test_tam_batch_traj_w1
tukey_table_tamanho_W1 <- tukey_table_tamanho_W1[, c("Source", setdiff(colnames(tukey_table_tamanho_W1), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_tamanho_w1.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Mean Size - w1", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_tamanho_W1)
dev.off()
#Tamanho na w2
tukey_table_tamanho_W2 <- tukey_test_tam_batch_traj_w2
tukey_table_tamanho_W2 <- tukey_table_tamanho_W2[, c("Source", setdiff(colnames(tukey_table_tamanho_W2), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_tamanho_w2.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Mean Size - w2", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_tamanho_W2)  # Desenha a tabela no arquivo PDF
dev.off()
#Porcentagem na w1
tukey_table_porcentagem_W1 <- tukey_test_porc_batch_traj_w1
tukey_table_porcentagem_W1 <- tukey_table_porcentagem_W1[, c("Source", setdiff(colnames(tukey_table_porcentagem_W1), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_porcentagem_w1.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Percentage of Small EV's - w1", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_porcentagem_W1)  # Desenha a tabela no arquivo PDF
dev.off()
#Porcentagem na w2
tukey_table_porcentagem_W2 <- tukey_test_porc_batch_traj_w2
tukey_table_porcentagem_W2 <- tukey_table_porcentagem_W2[, c("Source", setdiff(colnames(tukey_table_porcentagem_W2), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_porcentagem_w2.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Percentage of Small EV's - w2", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_porcentagem_W2)  # Desenha a tabela no arquivo PDF
dev.off()
#Concentração na w1
tukey_table_concentracao_W1 <- tukey_test_conc_batch_traj_w1
tukey_table_concentracao_W1 <- tukey_table_concentracao_W1[, c("Source", setdiff(colnames(tukey_table_concentracao_W1), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_concentracao_w1.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Concentration - w1", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_concentracao_W1)  # Desenha a tabela no arquivo PDF
dev.off()
#Concentração na w2
tukey_table_concentracao_W2 <- tukey_test_conc_batch_traj_w2
tukey_table_concentracao_W2 <- tukey_table_concentracao_W2[, c("Source", setdiff(colnames(tukey_table_concentracao_W2), "Source"))] #Source na 1a coluna
pdf("plots/tukey_table_concentracao_w2.pdf", width = 15, height = 25)
grid.newpage()  # Inicia uma nova página para o gráfico
grid.text("Tukey Test - Concentration - w2", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))  # Título no topo da página
grid.table(tukey_table_concentracao_W2)  # Desenha a tabela no arquivo PDF
dev.off()

## Teste de Dunn ajustado para batch
########Tamanho w1
for (b in unique(nanosight_sem_outliers_mean_w1$batch_nan)) {
  print(paste("Results for mean considering batch - w1", b))
  # Subset dos dados por batch
  subset_data_mean_w1 <- nanosight_sem_outliers_mean_w1[nanosight_sem_outliers_mean_w1$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_mean_w1 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_mean_w1)
  print(kw_result_mean_w1)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_mean_w1$p.value < 0.05) {
    dunn_result_mean_w1 <- dunn.test(subset_data_mean_w1$zscore_mean, g = subset_data_mean_w1$Trajetoria, method = "bonferroni")
    print(dunn_result_mean_w1)
  }
}
########Tamanho w2
for (b in unique(nanosight_sem_outliers_mean_w2$batch_nan)) {
  print(paste("Results for mean considering batch - w2", b))
  # Subset dos dados por batch
  subset_data_mean_w2 <- nanosight_sem_outliers_mean_w2[nanosight_sem_outliers_mean_w2$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_mean_w2 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_mean_w2)
  print(kw_result_mean_w2)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_mean_w2$p.value < 0.05) {
    dunn_result_mean_w2 <- dunn.test(subset_data_mean_w2$zscore_mean, g = subset_data_mean_w2$Trajetoria, method = "bonferroni")
    print(dunn_result_mean_w2)
  }
}
#Porcentagem w1
for (b in unique(nanosight_sem_outliers_porcentagem_w1$batch_nan)) {
  print(paste("Results for percentage considering batch - w1", b))
  # Subset dos dados por batch
  subset_data_porcentagem_w1 <- nanosight_sem_outliers_porcentagem_w1[nanosight_sem_outliers_porcentagem_w1$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_porcentagem_w1 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_porcentagem_w1)
  print(kw_result_porcentagem_w1)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_porcentagem_w1$p.value < 0.05) {
    dunn_result_porcentagem_w1 <- dunn.test(subset_data_porcentagem_w1$zscore_mean, g = subset_data_porcentagem_w1$Trajetoria, method = "bonferroni")
    print(dunn_result_porcentagem_w1)
  }
}
#Porcentagem w2
for (b in unique(nanosight_sem_outliers_porcentagem_w2$batch_nan)) {
  print(paste("Results for percentage considering batch - w2", b))
  # Subset dos dados por batch
  subset_data_porcentagem_w2 <- nanosight_sem_outliers_porcentagem_w2[nanosight_sem_outliers_porcentagem_w2$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_porcentagem_w2 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_porcentagem_w2)
  print(kw_result_porcentagem_w2)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_porcentagem_w2$p.value < 0.05) {
    dunn_result_porcentagem_w2 <- dunn.test(subset_data_porcentagem_w2$zscore_mean, g = subset_data_porcentagem_w2$Trajetoria, method = "bonferroni")
    print(dunn_result_porcentagem_w2)
  }
}
#Concentração w1
for (b in unique(nanosight_sem_outliers_concentracao_w1$batch_nan)) {
  print(paste("Results for concentration considering batch - w1", b))
  # Subset dos dados por batch
  subset_data_concentracao_w1 <- nanosight_sem_outliers_concentracao_w1[nanosight_sem_outliers_concentracao_w1$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_concentracao_w1 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_concentracao_w1)
  print(kw_result_concentracao_w1)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_concentracao_w1$p.value < 0.05) {
    dunn_result_concentracao_w1 <- dunn.test(subset_data_concentracao_w1$zscore_mean, g = subset_data_concentracao_w1$Trajetoria, method = "bonferroni")
    print(dunn_result_concentracao_w1)
  }
}
#Concentração w2
for (b in unique(nanosight_sem_outliers_concentracao_w2$batch_nan)) {
  print(paste("Results for concentration considering batch - w2", b))
  # Subset dos dados por batch
  subset_data_concentracao_w2 <- nanosight_sem_outliers_concentracao_w2[nanosight_sem_outliers_concentracao_w2$batch_nan == b, ]
  # Teste de Kruskal-Wallis
  kw_result_concentracao_w2 <- kruskal.test(zscore_mean ~ Trajetoria, data = subset_data_concentracao_w2)
  print(kw_result_concentracao_w2)
  # Teste de Dunn se o Kruskal-Wallis for significativo
  if (kw_result_concentracao_w2$p.value < 0.05) {
    dunn_result_concentracao_w2 <- dunn.test(subset_data_concentracao_w2$zscore_mean, g = subset_data_concentracao_w2$Trajetoria, method = "bonferroni")
    print(dunn_result_concentracao_w2)
  }
}
####Gerar tabelas de Kruskall-Wallis e Dunn (quando p < 0,05)
# Função para rodar Kruskal-Wallis e Dunn e salvar tabelas
run_tests_and_save <- function(data, file_prefix, variable_name) {
  for (b in unique(data$batch_nan)) {
    subset_data <- data[data$batch_nan == b, ]
    # Teste de Kruskal-Wallis
    kw_result <- kruskal.test(zscore_mean ~ Trajetoria, data)
    print(paste("Results for", variable_name, b))
    print(kw_result)
    # Salvar resultados de Kruskal-Wallis em uma tabela
    kw_table <- data.frame(
      Batch = b,
      Chi_Squared = kw_result$statistic,
      DF = kw_result$parameter,
      P_Value = kw_result$p.value
    )
    kw_table %>%
      gt() %>%
      gtsave(paste0("plots/", file_prefix, "_kruskal_", b, ".pdf"))
    # Teste de Dunn se o Kruskal-Wallis for significativo
    if (kw_result$p.value < 0.05) {
      dunn_result <- dunn.test(subset_data$zscore_mean, g = subset_data$Trajetoria, method = "bonferroni")
      print(dunn_result)
      # Transformar os resultados do Dunn em tabela
      dunn_table <- data.frame(
        Comparison = dunn_result$comparisons,
        Z = dunn_result$Z,
        `P-Value` = dunn_result$P,
        `Adjusted P-Value` = dunn_result$P.adjusted
      )
      # Exportar a tabela de Dunn Test
      dunn_table %>%
        gt() %>%
        gtsave(paste0("plots/", file_prefix, "_dunn_", b, ".pdf"))
    }
  }
}
     # Rodar a função para cada conjunto de dados e variável
run_tests_and_save(nanosight_sem_outliers_mean_w1, "mean_w1", "Mean W1")
run_tests_and_save(nanosight_sem_outliers_mean_w2, "mean_w2", "Mean W2")
run_tests_and_save(nanosight_sem_outliers_porcentagem_w1, "porcentagem_w1", "Percentage W1")
run_tests_and_save(nanosight_sem_outliers_porcentagem_w2, "porcentagem_w2", "Percentage W2")
run_tests_and_save(nanosight_sem_outliers_concentracao_w1, "concentracao_w1", "Concentration W1")
run_tests_and_save(nanosight_sem_outliers_concentracao_w2, "concentracao_w2", "Concentration W2")

#teste Pearson
cor.test(nanosight_plus_sampleinfo$concentracao_real, nanosight_plus_sampleinfo$tamanho_mean_average, method = "pearson") #0.27 (correlação media)
cor.test(nanosight_plus_sampleinfo_sem_outliers_mean_concentracao$tamanho_mean_average.x, nanosight_plus_sampleinfo_sem_outliers_mean_concentracao$concentracao_real.x, method = "pearson") #0.3363561

#teste de Spearman
cor.test(
  nanosight_plus_sampleinfo$concentracao_real, 
  nanosight_plus_sampleinfo$tamanho_mean_average, 
  method = "spearman"
) #0.433
cor.test(
  nanosight_plus_sampleinfo_sem_outliers_mean_concentracao$concentracao_real.x, 
  nanosight_plus_sampleinfo_sem_outliers_mean_concentracao$tamanho_mean_average.x, 
  method = "spearman"
) #0.422
