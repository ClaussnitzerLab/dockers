
devtools::install_github("cardiomoon/webr")

library(dplyr)
library(ggplot2)
library(lme4)
library(rlang)
library(qvalue)
library(readr)
library(moonBook)  # ?
library(webr)
library(ComplexHeatmap)
library(circlize)
library(vcfR)
library(ggsignif)
library(ggpubr)
library(stringr)
library(rstatix)  # 
rm(list = ls())

### plot gene expressions
plot_expressions <- function(plot_data, depot, output_path){
  figure_path_sex <- sprintf("%s/expression_day_sex_%s.pdf", output_path, depot)
  figure_path <- sprintf("%s/expression_day_%s.pdf", output_path, depot)
  plot_data["Day"] <- factor(plot_data$Day, levels = c(0, 3, 8, 14))
  Male_len <- dim(plot_data[plot_data$sex == 1,])[1]
  Female_len <- dim(plot_data[plot_data$sex == 2,])[1]
  plot_data["sex"] <- as.factor(ifelse(plot_data$sex==1, sprintf("Male (N=%d)", Male_len), sprintf("Female (N=%d)", Female_len)))
  plot_data["expression"] <- plot_data[gene_ensb_id]
  summary(plot_data$sex)
  summary(plot_data$Day)
  
  
  
  ggplot(plot_data, aes(x=Day, y=expression, fill=sex)) + geom_boxplot() + ggtitle(sprintf("Depot: %s", depot)) + xlab("Differentiation day") + ylab(sprintf("Expression %s(%s)", gene_name, gene_ensb_id))
  ggsave(figure_path_sex)
  
  ggplot(plot_data, aes(x=Day, y=expression)) +  geom_boxplot() + ggtitle(sprintf("Depot: %s (N=%d)", depot, dim(plot_data)[1])) + xlab("Differentiation day") + ylab(sprintf("Expression %s(%s)", gene_name, gene_ensb_id))
  ggsave(figure_path)
}

run_on_a_gene_for_plotting_expressions <- function(data, output_path, gene_name, gene_ensb_id){
  # create dir's
  # output_path <- sprintf("%s/%s", root_output_path, gene_name)
  dir.create(file.path(output_path))
  dir.create(file.path(output_path, 'gene_feature'))
  dir.create(file.path(output_path, 'heatmaps_qval'))
  dir.create(file.path(output_path, 'pie_chart_qval'))
  dir.create(file.path(output_path, 'eQTL'))
  
  #### get data
  curr_gene_expression_data <- select(data, c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", gene_ensb_id))
  
  #### remove FFA
  num_entries <- dim(curr_gene_expression_data)[1]
  curr_gene_expression_data <- curr_gene_expression_data[curr_gene_expression_data$FFA == 0,]
  writeLines(sprintf("%d enties removed due to FFA == 1 [current len: %d]", num_entries-dim(curr_gene_expression_data)[1], dim(curr_gene_expression_data)[1]))
  
  curr_vc_data <- curr_gene_expression_data[curr_gene_expression_data$cellType == "vc",]
  curr_sc_data <- curr_gene_expression_data[curr_gene_expression_data$cellType == "sc",]
  
  depot <- "Visceral"
  plot_data <- curr_vc_data
  plot_expressions(plot_data, depot, output_path)
  
  
  depot <- "Subcutaneous"
  plot_data <- curr_sc_data
  plot_expressions(plot_data, depot, output_path)
}

prepare_data_for_expression_imaging <- function(data, gene_ensb_id){
  #### regression section ####
  curr_lp_data <- data
  curr_lp_data["expression"] <- curr_lp_data[gene_ensb_id]
  
  #### remove FFA
  num_entries <- dim(curr_lp_data)[1]
  curr_lp_data <- curr_lp_data[curr_lp_data$FFA == 0,]
  writeLines(sprintf("%d enties removed due to FFA == 1 [current len: %d]", num_entries-dim(curr_lp_data)[1], dim(curr_lp_data)[1]))
  
  ### remove gene expressions except the queried one
  curr_lp_data <- curr_lp_data %>% select(-contains("ENSG"))
  
  ### get meta data 
  lp_meta <- select(curr_lp_data, c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))
  lp_imaging <- select(curr_lp_data, -c("X","batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))
  
  ### cleaning imaging features
  lp_imaging[,] <- sapply(sapply(lp_imaging[,], as.numeric), as.numeric)
  #######remove features by name you don't want#####
  lp_imaging<-as.data.frame(t(lp_imaging))
  lp_imaging<-lp_imaging[!grepl( "_Costes_" , rownames(lp_imaging) ) ,  ]
  lp_imaging<-lp_imaging[!grepl( "_Manders_" , rownames(lp_imaging)) ,  ]
  lp_imaging<-lp_imaging[!grepl( "_RWC_" , rownames(lp_imaging)) ,  ]
  lp_imaging<-lp_imaging[!grepl( "SmallBODIPY" , rownames(lp_imaging) ) ,  ]
  ######remove 0 features#######
  mean<-rowMeans(lp_imaging)
  var<-cbind(mean, lp_imaging)
  lp_imaging<-subset(var, var$mean != "0")
  ######remove NA features#####
  lp_imaging<- lp_imaging[complete.cases(lp_imaging), ]
  writeLines(sprintf("lp_imaging; [%d, %d ]", dim(lp_imaging)[1], dim(lp_imaging)[2]))
  lp_imaging<-as.data.frame(t(lp_imaging[,2:dim(lp_imaging)[2]]))
  curr_lp_data<-cbind(lp_meta, lp_imaging)
  curr_lp_data <- curr_lp_data[complete.cases(curr_lp_data), ]  
  return(curr_lp_data)
}

update_regression_arrays <- function(current_model_input, outcome, exposure, sex_label){
  # exposure <- "expression"
  # sex_label <- "Male"  "both_sex"
  # current_model_input <- current_model_input_male
  current_model_input$sex <- as.factor(current_model_input$sex)
  current_model_input$T2D <- as.factor(current_model_input$T2D)
  current_model_input$batch <- as.factor(current_model_input$batch)
  current_model_input$expression <- current_model_input$expression + .001
  # temp <- select(current_model_input, c("expression", "age", "BMI", "sex", "T2D", "batch"))
  # summary(temp)
  current_model_input <- na.omit(current_model_input)
  if (sex_label == "both_sex"){
    lm_model <- summary.lm(glm(get(outcome) ~ log(get(exposure))
                               + age
                               + BMI
                               +sex
                               +batch, 
                               data = current_model_input))
  }else{  # sex stratified, not adjusting for sex
    lm_model <- summary.lm(lm(get(outcome) ~ log(get(exposure))
                              + age
                              + BMI
                              +batch , 
                              data = current_model_input))
  }
  
  beta <-lm_model$coefficients[2, 1]
  pval <-lm_model$coefficients[2, 4]
  
  output <- list(beta, 0, pval, FALSE)
  return(output)
}

run_regression <- function(current_model_input, output_path, depot, day, sex_label, gene_ensb_id){
  imaging_data <- select(current_model_input, -c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))
  
  # sex_label <- "Male"
  outcome <- "Cells_Children_LargeBODIPYObjects_Count"
  output_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(output_df) <- c("variable", "coeff", "pvalue")
  for (outcome in colnames(imaging_data)){
    #writeLines(sprintf("%s", outcome))
    output <- update_regression_arrays(current_model_input, outcome, "expression", sex_label)
    out <- as.data.frame(outcome)
    colnames(out) <- c("variable")
    out["beta"] <- as.data.frame(unlist(output[1]))
    se <- as.array(unlist(output[2]))
    out["pval"] <- unlist(output[3])
    flag <- unlist(output[4])
    if (flag[1] == FALSE){
      output_df <- rbind(output_df, out)
    }
  }
  output_df <- na.omit(output_df)
  #########combine outcome and 
  p<-output_df$pval
  if (is_empty(p)){
    model_output_final <- output_df
    c <- c(colnames(model_output_final), "qvalues")
  }else{
    qvalues<-qvalue(p)
    q_value<-qvalues[["qvalues"]]
    
    model_output_final<-cbind(output_df, q_value)  
  }
  ############### SAVE
  write_csv(as.data.frame(model_output_final), sprintf("%s/%s_%d_%s.csv", output_path, depot, day, sex_label))  
}

running_regression_commands <- function(output_path, curr_lp_data, gene_ensb_id, day){
  ### run regression commands
  depots <- c("vc", "sc")
  sex_stratifications <- c(TRUE, FALSE) 
  gene_feature_output_path = sprintf("%s/gene_feature/", output_path)
  dir.create(gene_feature_output_path)
  for (depot in depots){
    for (sex_stratified in sex_stratifications){
      current_model_input<-subset(curr_lp_data, curr_lp_data$cellType == depot & curr_lp_data$Day == day)
      if (sex_stratified){
        current_model_input_male <- subset(current_model_input, current_model_input$sex == 1)
        current_model_input_female <- subset(current_model_input, current_model_input$sex == 2)
        run_regression(current_model_input_male, gene_feature_output_path, depot, day, "Male", gene_ensb_id)
        run_regression(current_model_input_female, gene_feature_output_path, depot, day, "Female", gene_ensb_id)
      }else{
        # no sex stratification, use the current_model_input
        run_regression(current_model_input, gene_feature_output_path, depot, day, "both_sex", gene_ensb_id)
      }
    }
  }

}

#### Pie chart ###
plot_pie_chart <- function(input_data, output_path, depot, day, sex_label){
  input_data["features"] <- input_data$variable
  ####compartment###
  pie_1 <- input_data %>% dplyr::mutate(Features = ifelse(grepl("Cells_", input_data$features),'Cells', 
                                                          ifelse(grepl("Cytoplasm_", input_data$features),'Cytoplasm', 
                                                                 ifelse(grepl("Nuclei", input_data$features),'Nuclei',
                                                                        'other_feature'))))
  ####channel###
  pie_2 <- pie_1 %>% dplyr::mutate(feature_org = ifelse(grepl("BODIPY", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'BODIPY', ifelse(grepl("BODIPY", pie_1$features) & grepl("AGP", pie_1$features), 'BODIPY_AGP', ifelse(grepl("BODIPY", pie_1$features) & grepl("Mito", pie_1$features), 'BODIPY_Mito', ifelse(grepl("BODIPY", pie_1$features) & grepl("DNA", pie_1$features),'BODIPY_DNA', ifelse(grepl("Mito", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("DNA", pie_1$features), 'Mito', ifelse(grepl("Mito", pie_1$features)  & grepl("AGP", pie_1$features) ,'Mito_AGP',  ifelse(grepl("Mito", pie_1$features)  & grepl("DNA", pie_1$features) ,'Mito_DNA', ifelse(grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'AGP', ifelse(grepl("AGP", pie_1$features)  & grepl("DNA", pie_1$features) ,'AGP_DNA', ifelse(grepl("DNA", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("BODIPY", pie_1$features), 'DNA', 'other_feature')))))))))))
  
  ####measurement###
  pie_3 <- pie_2 %>% dplyr::mutate(feature_measure = 
                                     ifelse(grepl("Texture", pie_2$features),'Texture', 
                                            ifelse(grepl("Intensity", pie_2$features),'Intensity', 
                                                   ifelse(grepl("Granularity", pie_2$features),'Granularity',
                                                          'other_feature'))))
  
  ###plot pie charts###
  
  
  ###plot compartments and channels
  pdf(sprintf("%s/%s_%d_%s.pdf", output_path, depot, day, sex_label))  # 
  PieDonut(pie_3,aes(pies= Features , donuts =  feature_org), title=sprintf("%s-Day:%d;sex:%s [N=%d]", depot, day, sex_label, dim(input_data)[1]))
  dev.off()
}

try_plot <- function(curr_plot_data, output_path, depot, day, sex_label){
  curr_plot_data <- subset(curr_plot_data, curr_plot_data$q_value < .05)
  if (dim(curr_plot_data)[1] > 10){
    plot_pie_chart(curr_plot_data, output_path, depot, day, sex_label)
  }
}

generate_pie_chart <- function(output_path, day){

  depots <- c("vc", "sc")
  gene_feature_output_path = sprintf("%s/gene_feature/", output_path)
  gene_feature_piechart_output_path = sprintf("%s/pie_chart_qval/", output_path)
  dir.create(gene_feature_piechart_output_path)
  for (depot in depots){
    Male_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_output_path, depot, day, "Male"))
    try_plot(Male_data, gene_feature_piechart_output_path, depot, day, "Male")
    
    Female_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_output_path, depot, day, "Female"))
    try_plot(Female_data, gene_feature_piechart_output_path, depot, day, "Female")
    
    both_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_output_path, depot, day, "both_sex"))
    try_plot(both_data, gene_feature_piechart_output_path, depot, day, "both_sex")
  }

}  

#### heatmap plots ####
plot_heatmap <- function(nonredundant_significant_features, output_image_path){
  pie_1 <- nonredundant_significant_features %>% dplyr::mutate(feature_location = ifelse(grepl("Cells_", nonredundant_significant_features$features),'Cells', ifelse(grepl("Cytoplasm_", nonredundant_significant_features$features),'Cytoplasm', ifelse(grepl("Nuclei", nonredundant_significant_features$features),'Nuclei', 'other_feature'))))
  
  ####channel###
  pie_2 <- pie_1 %>% dplyr::mutate(feature_org =
                                     ifelse(grepl("BODIPY", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'BODIPY',
                                            ifelse(grepl("BODIPY", pie_1$features) & grepl("AGP", pie_1$features), 'BODIPY_AGP',
                                                   ifelse(grepl("BODIPY", pie_1$features) & grepl("Mito", pie_1$features), 'BODIPY_Mito',
                                                          ifelse(grepl("BODIPY", pie_1$features) & grepl("DNA", pie_1$features),'BODIPY_DNA',
                                                                 ifelse(grepl("Mito", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("DNA", pie_1$features), 'Mito',
                                                                        ifelse(grepl("Mito", pie_1$features)  & grepl("AGP", pie_1$features) ,'Mito_AGP',
                                                                               ifelse(grepl("Mito", pie_1$features)  & grepl("DNA", pie_1$features) ,'Mito_DNA',
                                                                                      ifelse(grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'AGP',
                                                                                             ifelse(grepl("AGP", pie_1$features)  & grepl("DNA", pie_1$features) ,'AGP_DNA',
                                                                                                    ifelse(grepl("DNA", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("BODIPY", pie_1$features), 'DNA',
                                                                                                           'other_feature')))))))))))
  
  ####measurement###
  pie_3 <- pie_2 %>% dplyr::mutate(feature_measure =  ifelse(grepl("Texture", pie_2$features),'Texture',  ifelse(grepl("Intensity", pie_2$features),'Intensity', ifelse(grepl("Granularity", pie_2$features),'Granularity', 'other_feature'))))
  
  pie_4 <- pie_3 %>% dplyr::mutate(significant = ifelse(pie_3$q_value <= 0.05 & pie_3$beta < 0,'negative', ifelse(pie_3$q_value <= 0.05 & pie_3$beta > 0,'positive', 'other_feature')))
  
  matrix<-as.matrix(pie_4$beta)
  f2 <- colorRamp2(seq(min(matrix), max(matrix), length = 6), c( "#fafafa","#f2d38a","#f2bf80", "#f29180","#f29180","#f01d1d"))
  pdf(output_image_path)
  plot = Heatmap(matrix, col = f2, name = "effect size", show_row_names = F,  column_names_gp = gpar(col = "black", fontsize = 2),
                 right_annotation = rowAnnotation(channel =  pie_4$feature_org, measurment = pie_4$feature_measure, significance = pie_4$significant,
                                                  col= list(channel = c("Mito" = "#f56464" , "AGP" = "#f2cf41", "BODIPY" = "#4dac26", "DNA" = "#3b8fe3", "other_feature" = "#959ca3",
                                                                        "AGP_DNA" = "#d7dce0", "BODIPY_AGP" = "#d7dce0", "Mito_AGP" = "#d7dce0", "BODIPY_Mito" = "#d7dce0", "BODIPY_DNA" = "#d7dce0",
                                                                        "Mito_DNA" = "#d7dce0"), measurment =c("Intensity" = "#b5bac7" , "Granularity" = "#626a80", "Texture" = "#343f5c", 
                                                                                                               "other_feature" = "#e4e7f0"),  significance = c("negative"= "#9BC7E4",  "positive"= "#B02636", 
                                                                                                                                                               "other_feature" = "white"))))
  
  draw(plot)
  dev.off()
  
}

try_heatmap <- function(curr_data, curr_data_output_path){
  curr_data <- subset(curr_data, curr_data$q_value < .05)
  curr_data["features"] <- curr_data$variable
  if (dim(curr_data)[1] > 10){
    plot_heatmap(curr_data, curr_data_output_path)
    
  }  
}

generate_heatmap_plots <- function(output_path, day){
  depots <- c("vc", "sc")
  gene_feature_input_path = sprintf("%s/gene_feature/", output_path)
  gene_feature_output_path = sprintf("%s/heatmaps_qval/", output_path)
  dir.create(gene_feature_output_path)
  depot <- "sc"
  for (depot in depots){
    # Male
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "Male"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "Male")
    try_heatmap(curr_data, curr_data_output_path)
    
    
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "Female"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "Female")
    try_heatmap(curr_data, curr_data_output_path)
    
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "both_sex"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "both_sex")
    try_heatmap(curr_data, curr_data_output_path)
  }

}


try_volcano <- function(curr_data, curr_data_output_path, volcano_title){
  if (dim(curr_data)[1] < 10){
    return()
  }
  th <- .05
  lower_fdr_variable <- sprintf("qvalue>=%.02f", th)
  higher_fdr_variable <- sprintf("qvalue<%.02f", th)
  signif_variable <- sprintf("qvalue<%.02f, beta>.5", th)
  new_colors <- c("blue", "blue4", "red")
  names(new_colors) <- c(lower_fdr_variable, higher_fdr_variable, signif_variable)
  
  curr_data$ColorCode <- lower_fdr_variable  # every one gets low fdr color
  curr_data$ColorCode[curr_data$q_value <= th] <- higher_fdr_variable  # the onese with good fdr get a color
  curr_data$ColorCode[abs(curr_data$beta) > 0.5 & curr_data$q_value <= th] <- signif_variable  # significant ones
 
  p <- ggplot(data=curr_data, aes(x=beta, y=-log10(q_value), col=ColorCode)) + geom_point() +
    theme_minimal() + 
    scale_colour_manual(values = new_colors)+
    labs(title=sprintf("%s", volcano_title), 
         x ="beta", y = "-log10(qvalue)") + 
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text = element_text(size = 18),
      # legend.position = "none",
      plot.background = element_rect(fill = "white")
    )
  p
  ggsave(curr_data_output_path)
}

generate_volcano_chart <- function(output_path, day){
  depots <- c("vc", "sc")
  gene_feature_input_path = sprintf("%s/gene_feature/", output_path)
  gene_feature_output_path = sprintf("%s/volcano_qval/", output_path)
  dir.create(gene_feature_output_path)
  depot <- "sc"
  for (depot in depots){
    # Male
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "Male"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "Male")
    volcano_title <- sprintf("features %s;%s;%s", depot, day, "Male")
    try_volcano(curr_data, curr_data_output_path, volcano_title)
    
    
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "Female"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "Female")
    volcano_title <- sprintf("features %s;%s;%s", depot, day, "Female")
    try_volcano(curr_data, curr_data_output_path, volcano_title)
    
    curr_data <- read.csv(sprintf("%s/%s_%d_%s.csv", gene_feature_input_path, depot, day, "both_sex"))
    curr_data_output_path <- sprintf("%s/%s_%d_%s.pdf", gene_feature_output_path, depot, day, "both_sex")
    volcano_title <- sprintf("features %s;%s;%s", depot, day, "both_sex")
    try_volcano(curr_data, curr_data_output_path, volcano_title)
  }  
}


### For a gene ###
setwd(getwd())
args = commandArgs(trailingOnly=TRUE)
LP_data_path <- args[1]
gene_ensb_id <- args[2]
gene_name <- args[3]
day <- as.numeric(args[4])
output_path <- "LP_output"
dir.create(output_path)

data<-read.csv(file=LP_data_path)

# plotting expressions
if (day == 0){
  run_on_a_gene_for_plotting_expressions(data, output_path, gene_name, gene_ensb_id)  
}

# running regression models
curr_lp_data <- prepare_data_for_expression_imaging(data, gene_ensb_id)
running_regression_commands(output_path, curr_lp_data, gene_ensb_id, day)

# generating heatmap
generate_heatmap_plots(output_path, day)

# genrate piechart
generate_pie_chart(output_path, day)

# generate volcano plot
generate_volcano_chart(output_path, day)
