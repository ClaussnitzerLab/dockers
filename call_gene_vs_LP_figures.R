
library(dplyr)
library(ggplot2)
library(lme4)
library(rlang)
library(qvalue)
library(readr)
# library(ComplexHeatmap)
library(circlize)
library(vcfR)
library(ggsignif)
library(ggpubr)
library(stringr)
library(tidyr)
library(cowplot)
rm(list = ls())

# For reporting stats
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = sprintf("N=%d", length(y))))
}

#### plotting the expression level of a gene across the differnetiation time points ####
plot_expressions <- function(plot_data, depot, ffa_tag, output_path){
  figure_path_sex <- sprintf("%s/expression_day_sex_%s_%s.pdf", output_path, depot, ffa_tag)
  figure_path <- sprintf("%s/expression_day_%s_%s.pdf", output_path, depot, ffa_tag)
  csv_path <- sprintf("%s/expression_day_%s_%s.csv", output_path, depot, ffa_tag)
  plot_data["Day"] <- factor(plot_data$Day, levels = c(0, 3, 8, 14))
  
  ll <- plot_data[plot_data$sex == 1,]
  Male_len <- dim(plot_data[plot_data$sex == 1,])[1]
  Female_len <- dim(plot_data[plot_data$sex == 2,])[1]
  # plot_data["sex"] <- as.factor(ifelse(plot_data$sex==1, sprintf("Male (N=%d)", Male_len), sprintf("Female (N=%d)", Female_len)))
  plot_data["sex"] <- as.factor(ifelse(plot_data$sex==1, sprintf("Male"), sprintf("Female")))
  plot_data["expression"] <- plot_data[gene_ensb_id]
  summary(plot_data$sex)
  summary(plot_data$Day)
  write.csv(plot_data, csv_path, row.names = FALSE)
  ggplot(plot_data, aes(x=Day, y=expression, fill=sex)) + geom_boxplot() + 
    ggtitle(sprintf("Depot: %s", depot)) + 
    xlab("Differentiation day") + ylab(sprintf("Expression %s(%s)", gene_name, gene_ensb_id))+ 
    theme(text = element_text(size = 28)) + theme(legend.position = c(0.8, 0.8)) + 
    stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, position=position_dodge(.9)) +
    theme_classic()
  ggsave(figure_path_sex)
  
  ggplot(plot_data, aes(x=Day, y=expression)) +  geom_boxplot() + 
    ggtitle(sprintf("Depot: %s (N=%d)", depot, dim(plot_data)[1])) + 
    xlab("Differentiation day") + ylab(sprintf("Expression %s(%s)", gene_name, gene_ensb_id))+ 
    theme(text = element_text(size = 20))   + theme(legend.position = c(0.8, 0.8))  
  ggsave(figure_path)
}

# preparing data for plotting gene expressions across differentiation time points
plotting_expression_vs_timepoints <- function(data, output_path, gene_name, gene_ensb_id){
  curr_gene_expression_data <- select(data, c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", gene_ensb_id))
  for (cellType in c("vc", "sc")){
    curr_cellType_data <- curr_gene_expression_data[curr_gene_expression_data$cellType == cellType,]
    if (cellType == "vc"){
      depot <- "Visceral"
    }else{
      depot <- "Subcutaneous"
    }
    for (ffa in c(0, 1)){
      if (ffa == 0){
        ffa_tag <- "FFA0"
      }else{
        ffa_tag <- "FFA1"
      }
      curr_ffa_data <- curr_cellType_data[curr_cellType_data$FFA == ffa,]
      plot_expressions(curr_ffa_data, depot, ffa_tag, output_path)
    }    
  }
}

# a sub-function for the regression function. This generates a row of the output table
update_regression_arrays <- function(current_model_input, outcome, exposure, sex_label){
  # exposure <- "expression"
  # sex_label <- "Male"  "both_sex"
  # current_model_input <- current_model_input_male
  current_model_input$sex <- as.factor(current_model_input$sex)
  current_model_input$T2D <- as.factor(current_model_input$T2D)
  current_model_input$batch <- as.factor(current_model_input$batch)
  current_model_input$expression <- current_model_input$expression + .001
  current_model_input <- na.omit(current_model_input)

  formula <- sprintf("get('%s') ~ log(get('%s')) + age", outcome, exposure)
  if (sex_label == "both_sex"){
    if (length(unique(current_model_input$sex)) > 1){
      formula <- sprintf("%s + sex", formula)
    }        
  }
  if (length(unique(current_model_input$BMI)) > 1){
    formula <- sprintf("%s + BMI", formula)
  }
  #if (length(unique(current_model_input$T2D)) > 1){
  #  formula <- sprintf("%s + T2D", formula)
  #}
  if (length(unique(current_model_input$batch)) > 1){
    formula <- sprintf("%s + batch", formula)
  } 
  lm_model <- summary.lm(lm(formula = formula, data = current_model_input))

  beta <-lm_model$coefficients[2, 1]
  pval <-lm_model$coefficients[2, 4]
  
  output <- list(beta, 0, pval, FALSE)
  return(output)
}



# running the lm function
run_regression <- function(current_model_input, output_path, depot, day, sex_label, FFA_tag, gene_ensb_id){
  imaging_data <- select(current_model_input, -c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))
  output_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(output_df) <- c("variable", "coeff", "pvalue")
  for (outcome in colnames(imaging_data)){
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
    if (length(p) > 20){
      qvalues<-qvalue(p)
      q_value<-qvalues[["qvalues"]]
    }else{
      qvalues<-qvalue(p, lambda=0)
      q_value<-qvalues[["qvalues"]]
    }
    
    model_output_final<-cbind(output_df, q_value)  
  }
  write_csv(as.data.frame(model_output_final), sprintf("%s/%s_%d_%s_%s.csv", output_path, depot, day, sex_label, FFA_tag))  
}

####### DONUT PLOT #######
add_compartment_names <- function(results_df){
  results_df <- results_df %>%
    dplyr::mutate(NCC_type = 
                    ifelse(grepl("Nuclei", rownames(results_df)) & !grepl("Cells|Cytoplasm)", rownames(results_df)),'Nuclei', 
                           ifelse(grepl("Cytoplasm", rownames(results_df)) & !grepl("Cells|Nuclei", rownames(results_df)),'Cytoplasm', 
                                  'Cells')))
  return(results_df)
}

add_channel_names <- function(results_df){
  results_df <- results_df %>% dplyr::mutate(feature_color = 
                                               ifelse(grepl("BODIPY", rownames(results_df)) & !grepl("AGP|Mito|DNA", rownames(results_df)),'Lipid',
                                                      ifelse(grepl("AGP", rownames(results_df)) & !grepl("BODIPY|Mito|DNA", rownames(results_df)),'AGP',
                                                             ifelse(grepl("Mito", rownames(results_df)) & !grepl("AGP|BODIPY|DNA", rownames(results_df)),'Mito',
                                                                    ifelse(grepl("DNA", rownames(results_df)) & !grepl("AGP|Mito|BODIPY", rownames(results_df)),'DNA',
                                                                           ifelse(!grepl("AGP|Mito|DNA|BODIPY", rownames(results_df)), 'Compartmental',
                                                                                  'Combination'))))))
  return(results_df)
}

add_combination_channel <- function(results_df){
  results_df <- results_df %>%
    dplyr::mutate(feature_color2 = 
                    ifelse(grepl("BODIPY", rownames(results_df)) & !grepl("AGP|Mito|DNA", rownames(results_df)),'Lipid',
                           ifelse(grepl("AGP", rownames(results_df)) & !grepl("BODIPY|Mito|DNA", rownames(results_df)),'AGP',
                                  ifelse(grepl("Mito", rownames(results_df)) & !grepl("AGP|BODIPY|DNA", rownames(results_df)),'Mito',
                                         ifelse(grepl("DNA", rownames(results_df)) & !grepl("AGP|Mito|BODIPY", rownames(results_df)),'DNA',
                                                ifelse(grepl("(?=.*BODIPY)(?=.*AGP)", rownames(results_df), perl = TRUE), 'Lipid/AGP',
                                                       ifelse(grepl("(?=.*BODIPY)(?=.*DNA)", rownames(results_df), perl = TRUE), 'Lipid/DNA',
                                                              ifelse(grepl("(?=.*BODIPY)(?=.*Mito)", rownames(results_df), perl = TRUE), 'Lipid/Mito',
                                                                     ifelse(grepl("(?=.*AGP)(?=.*Mito)", rownames(results_df), perl = TRUE), 'AGP/Mito',
                                                                            ifelse(grepl("(?=.*AGP)(?=.*DNA)", rownames(results_df), perl = TRUE), 'AGP/DNA',
                                                                                   ifelse(grepl("(?=.*DNA)(?=.*Mito)", rownames(results_df), perl = TRUE), 'DNA/Mito',
                                                                                          'Compartmental')))))))))))
  return(results_df)
}

add_measurement_names <- function(results_df){
  results_df <- results_df %>%
    dplyr::mutate(feature_type = 
                    ifelse(grepl("Intensity", rownames(results_df)) & !grepl("Texture|Granularity|Correlation|AreaShape|RadialDistribution|Neighbors|Count|Number", rownames(results_df)),'Intensity', #Counting Location_CenterMassIntensity and Location_MaxIntensity features as just intensity features
                           ifelse(grepl("Texture", rownames(results_df)) & !grepl("Intensity|Granularity|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", rownames(results_df)),'Texture', #Counting Texture_Correlation as just texture features
                                  ifelse(grepl("Granularity", rownames(results_df)) & !grepl("Texture|Intensity|Correlation|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", rownames(results_df)),'Granularity',
                                         ifelse(grepl("Correlation", rownames(results_df)) & !grepl("Texture|Granularity|Intensity|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", rownames(results_df)), 'Correlation',
                                                ifelse(grepl("RadialDistribution|Location|Neighbors", rownames(results_df)) & !grepl("Texture|Granularity|Correlation|Intensity|AreaShape|Count|Number", rownames(results_df)), 'Position', 'Shape/Size/Count')))))) #Counting "Cells_Parent_Nuclei" feature as Shape/Size/Count    
}

MyPieDonut <- function(results_df, in_title, file_name) {
  rownames(results_df) <- results_df$variable
  
  #Add a column with the compartment as a factor
  results_df <- add_compartment_names(results_df)
  results_df$NCC_type <- factor(results_df$NCC_type, levels = c("Nuclei", "Cytoplasm", "Cells"))
  
  #Add a column with the channel as a factor
  results_df <- add_channel_names(results_df)
  results_df$feature_color<-factor(results_df$feature_color, levels = c("Mito","AGP","Lipid","DNA", "Compartmental", "Combination"))
  
  
  #Add another column with the channel as a factor (specifying which combination for Combination features)
  results_df <- add_combination_channel(results_df)
  results_df$feature_color2<-factor(results_df$feature_color2, levels = c("Mito","AGP","Lipid","DNA", "Lipid/AGP", "Lipid/DNA", "Lipid/Mito", "AGP/Mito", "AGP/DNA", "DNA/Mito", "Compartmental"))
  
  #Add a column with the measurement as a factor
  #MEASUREMENTS: Intensity, Texture, Granularity, Position(Radial Distribution, Location, Neighbors), Shape/Size/Count(AreaShape, Count, Number), Correlation
  results_df <- add_measurement_names(results_df)
  results_df$feature_type <- factor(results_df$feature_type, levels = c("Correlation", "Intensity", "Texture", "Granularity", "Shape/Size/Count", "Position"))
  
  ##### start of Maya's scripta ####
  #####These are all of the categories for the features (color and type).
  color <- c("Mito","AGP","Lipid","DNA", "Combination", "Compartmental")
  type <- c("Intensity", "Texture", "Granularity", "Correlation", "Shape/Size/Count", "Position")
  
  #####slices_big is the dataframe that will be used to make the outer pie plot.
  #Make a table counting the number of features in each category.
  slices_big <- results_df %>% group_by(feature_color, feature_type) %>% summarise(feature_number = n())
  
  #Change all NA to zero.
  slices_big[is.na(slices_big)] <- 0
  #Name the columns of the data frame.
  colnames(slices_big) <- c("color", "type", "value")
  #Add a column that lists the color and type combination for each row.
  slices_big$category <- sprintf("%s_%s", slices_big$color, slices_big$type)
  #Specify the order that the slices of the pie charts will be plotted. This ensures that all of the types for a given color (Mito, AGP, etc.) are grouped together and align with the other pie chart.
  slices_big$category <- factor(slices_big$category, levels = rev(c("Mito_Correlation", "Mito_Intensity", "Mito_Texture", "Mito_Granularity", "Mito_Position", "Mito_Shape/Size/Count", "AGP_Correlation", "AGP_Intensity", "AGP_Texture", "AGP_Granularity", "AGP_Position", "AGP_Shape/Size/Count", "Lipid_Correlation", "Lipid_Intensity", "Lipid_Texture", "Lipid_Granularity", "Lipid_Position", "Lipid_Shape/Size/Count", "DNA_Correlation", "DNA_Intensity", "DNA_Texture", "DNA_Granularity", "DNA_Position", "DNA_Shape/Size/Count", "Combination_Correlation", "Combination_Intensity", "Combination_Texture", "Combination_Granularity", "Combination_Position", "Combination_Shape/Size/Count", "Compartmental_Correlation", "Compartmental_Intensity", "Compartmental_Texture", "Compartmental_Granularity", "Compartmental_Position", "Compartmental_Shape/Size/Count")))
  #Add a column with the percent of features in each category.
  slices_big$percent <- round(slices_big$value/sum(slices_big$value)*100, 1)
  #Add a column that specifies the labels for each slice of the pie.
  slices_big$labels <- paste0(as.character(slices_big$type), "\nn=", as.character(slices_big$value))
  #Remove the label for slices under a certain percent (2%). This will help to prevent labels from overcrowding and overlapping.
  for(u in 1:nrow(slices_big)) {
    if(slices_big[u,"percent"]<2.15) {
      slices_big[u,"labels"] <- ""
    }
  }
  #Reorder the rows of the dataframe so that the labels are centered in the correct slice. Without this, the labels may not be positioned correctly depending on the relative size of the slices.
  slices_big <- slices_big[match(c("Mito_Correlation", "Mito_Intensity", "Mito_Texture", "Mito_Granularity", "Mito_Position", "Mito_Shape/Size/Count", "AGP_Correlation", "AGP_Intensity", "AGP_Texture", "AGP_Granularity", "AGP_Position", "AGP_Shape/Size/Count", "Lipid_Correlation", "Lipid_Intensity", "Lipid_Texture", "Lipid_Granularity", "Lipid_Position", "Lipid_Shape/Size/Count", "DNA_Correlation", "DNA_Intensity", "DNA_Texture", "DNA_Granularity", "DNA_Position", "DNA_Shape/Size/Count", "Combination_Correlation", "Combination_Intensity", "Combination_Texture", "Combination_Granularity", "Combination_Position", "Combination_Shape/Size/Count", "Compartmental_Correlation", "Compartmental_Intensity", "Compartmental_Texture", "Compartmental_Granularity", "Compartmental_Position", "Compartmental_Shape/Size/Count"), slices_big$category),]
  
  #####plot_big is the plot of the outer pie chart.
  plot_big <- ggplot(slices_big, aes(x=2, y = value)) +
    geom_bar(stat="identity", width=1, color = "white", aes(fill = category)) +
    #This specifies the shades I chose. Intensity is 4 shades darker than the inner pie chart. Texture is 2 shades darker. Granularity is 2 shades lighter. Other is 4 shades lighter.
    scale_fill_manual(values = c("Mito_Correlation" = "#550505", "Mito_Intensity" = "#9e0a0a", "Mito_Texture" = "#e80f0f", "Mito_Granularity" = "#f45252", "Mito_Position" = "#f99b9b", "Mito_Shape/Size/Count" = "#fde5e5", "AGP_Correlation" = "#554606", "AGP_Intensity" = "#9f820b", "AGP_Texture" = "#e8bd10", "AGP_Granularity" = "#f3d453", "AGP_Position" = "#f8e69d", "AGP_Shape/Size/Count" = "#fdf9e6", "Lipid_Correlation" = "#224c11", "Lipid_Intensity" = "#3f8c1f", "Lipid_Texture" = "#5bcc2d", "Lipid_Granularity" = "#8cde6a", "Lipid_Position" = "#bdecaa", "Lipid_Shape/Size/Count" = "#effaea", "DNA_Correlation" = "#0c2335", "DNA_Intensity" = "#1b4c74", "DNA_Texture" = "#2a76b4", "DNA_Granularity" = "#559dd7", "DNA_Position" = "#a5cbea", "DNA_Shape/Size/Count" = "#e4f0f9", "Combination_Correlation" = "#2b1641", "Combination_Intensity" = "#522a7c", "Combination_Texture" = "#793db7", "Combination_Granularity" = "#a172d0", "Combination_Position" = "#c8ade4", "Combination_Shape/Size/Count" = "#f0e8f7", "Compartmental_Correlation" = "#2d3034", "Compartmental_Intensity" = "#51575e", "Compartmental_Texture" = "#767f88", "Compartmental_Granularity" = "#a0a6ac", "Compartmental_Position" = "#cacdd1", "Compartmental_Shape/Size/Count" = "#f4f4f5")) +
    #Position the labels in the center of each slice, outside of the pie chart.
    geom_text(aes(x=2.75, color = category, label = labels), size = 2.15, position = position_stack(vjust = 0.5)) +
    #Change the color of each label to match the fill color of its respective slice.
    scale_color_manual(values = c("Mito_Correlation" = "#550505", "Mito_Intensity" = "#9e0a0a", "Mito_Texture" = "#e80f0f", "Mito_Granularity" = "#f45252", "Mito_Position" = "#f99b9b", "Mito_Shape/Size/Count" = "#550505", "AGP_Correlation" = "#554606", "AGP_Intensity" = "#9f820b", "AGP_Texture" = "#e8bd10", "AGP_Granularity" = "#f3d453", "AGP_Position" = "#f8e69d", "AGP_Shape/Size/Count" = "#554606", "Lipid_Correlation" = "#224c11", "Lipid_Intensity" = "#3f8c1f", "Lipid_Texture" = "#5bcc2d", "Lipid_Granularity" = "#8cde6a", "Lipid_Position" = "#bdecaa", "Lipid_Shape/Size/Count" = "#224c11", "DNA_Correlation" = "#0c2335", "DNA_Intensity" = "#1b4c74", "DNA_Texture" = "#2a76b4", "DNA_Granularity" = "#559dd7", "DNA_Position" = "#a5cbea", "DNA_Shape/Size/Count" = "#0c2335", "Combination_Correlation" = "#2b1641", "Combination_Intensity" = "#522a7c", "Combination_Texture" = "#793db7", "Combination_Granularity" = "#a172d0", "Combination_Position" = "#c8ade4", "Combination_Shape/Size/Count" = "#2b1641", "Compartmental_Correlation" = "#2d3034", "Compartmental_Intensity" = "#51575e", "Compartmental_Texture" = "#767f88", "Compartmental_Granularity" = "#a0a6ac", "Compartmental_Position" = "#cacdd1", "Compartmental_Shape/Size/Count" = "#2d3034")) +
    #Turn it from a stacked bar plot into a pie plot.
    coord_polar("y", start=0) + 
    #Remove unnecessary formatting. Use blank theme.
    theme_void() +
    #Remove the legend
    theme(legend.position="none") +
    xlim(0.5, 2.75)
  
  #####slices_small is the dataframe that will be used to make the inner pie plot.
  #Make a table counting the number of features for each color (Mito, AGP, Lipid, DNA, Other).
  slices_small <- results_df %>% group_by(feature_color) %>% summarise(feature_number = n())
  #Change all NA to zero.
  slices_small[is.na(slices_small)] <- 0
  #Name the columns of the data frame.
  colnames(slices_small) <- c("color", "value")
  #Add a column with the percent of features for each color.
  slices_small$percent <- round(slices_small$value/sum(slices_small$value)*100, 1)
  #Specify the order that the slices of the pie charts will be plotted. This ensures that the two pie charts will align with one another.
  slices_small$color <- factor(slices_small$color, levels = rev(color))
  #Add a column that specifies the labels for each slice of the pie.
  slices_small$labels <- paste0(slices_small$color,"\n", "(", slices_small$percent, "%,", "\n", "n=", slices_small$value, ")" )
  #Remove the label for slices under a certain percent (3.5%). This will help to prevent labels from overcrowding and overlapping. This may not be necessary for the inner pie plot.
  for(w in 1:nrow(slices_small)) {
    if(slices_small[w,"percent"]<3.5) {
      slices_small[w,"labels"] <- ""
    }
  }
  #Garantees missing dyes do not compromise 
  color <- color[color %in% slices_small$color]
  #Reorder the rows of the dataframe so that the labels are centered in the correct slice. Without this, the labels may not be positioned correctly depending on the relative size of the slices.
  slices_small <- slices_small[match(color, slices_small$color),]
  
  #####plot_small is the plot of the inner pie chart.
  plot_small <- ggplot(slices_small, aes(x=2, y = value)) +
    geom_bar(stat="identity", width=1, color = "white", aes(fill = color)) + 
    #Specify the colors for Lipocyte Profiler. These are the hex codes that Sophie used.
    scale_fill_manual(values = c("Mito"="#f56464","AGP"="#f2cf41","Lipid"="#4dac26","DNA"="#65a6db", "Combination" = "#8d55c6", "Compartmental"="#959ca3")) +
    #Position the labels in the center of each slices towards the outer edge of the pie plot.
    geom_text(aes(x=2, label = labels), size = 2.5, position = position_stack(vjust = 0.5)) +
    #Turn it from a stacked bar plot into a pie plot.
    coord_polar("y", start=0) + 
    #Remove unnecessary formatting. Use blank theme.
    theme_void() +
    #Remove the legend
    theme(legend.position="none") +
    #Change the lower limit to adjust the size of the hole in the center of the plot.
    xlim(0.5, 2.5)
  
  #####Combine the plots.
  #Specify the size and position of the individual plots.
  size_big <- 1
  position_big <- 0
  size_small <- .6*size_big
  position_small <- position_big+size_big/2-size_small/2
  position_x_label <- 0.5
  position_y_label <- 0.95
  position_label2 <- 0.5*size_big
  #Overlay the pie plots and add titles.
  plot_combined <- ggdraw() +
    draw_plot(plot_big, x = position_big, y = position_big, width = size_big, height = size_big) +
    draw_plot(plot_small, x = position_small, y = position_small, width = size_small, height = size_small) +
    #Title at the top of the plot. Based on the file_name input to the function.
    draw_label(sprintf("%s", in_title), x = position_x_label, y = position_y_label, size = 12) +
    #Label in the hollow circle in the center of the plot.
    draw_label(paste0("Feature\nn = ", as.character(sum(slices_small$value))), x = position_label2, y = position_label2, size = 12)
  
  #####Export the final plot as a PDF.
  #Make sure that the root output directory has already been properly defined outside of the function.
  pdf(sprintf("%s.pdf", file_name))
  print(plot_combined)
  dev.off()
  
  write.csv(results_df, sprintf("%s.csv", file_name))
  return(plot_combined)
}
####### END of DONUT PLOT #######

### For a gene ###
setwd(getwd())
args = commandArgs(trailingOnly=TRUE)
LP_data_path <- args[1]
gene_ensb_id <- args[2]
gene_name <- args[3]
output_path <- "LP_output"
dir.create(output_path)

data<-read.csv(file=LP_data_path)
gene_feature_output_path = sprintf("%s/gene_feature/", output_path)
dir.create(gene_feature_output_path)

##########################
#### PLOT TIME POINTS ####
plotting_expression_vs_timepoints(data, output_path, gene_name, gene_ensb_id)

########################
#### RUN REGRESSION ####
curr_lp_data <- data
curr_lp_data["expression"] <- curr_lp_data[gene_ensb_id]
### remove gene expressions except the queried one
curr_lp_data <- curr_lp_data %>% select(-contains("ENSG"))

### get meta data 
lp_meta <- select(curr_lp_data, c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))
lp_imaging <- select(curr_lp_data, -c("X","batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression"))

### cleaning imaging features
lp_imaging[,] <- sapply(sapply(lp_imaging[,], as.numeric), as.numeric)

lp_imaging<-as.data.frame(t(lp_imaging))
lp_imaging<-lp_imaging[!grepl( "_Costes_" , rownames(lp_imaging) ) ,  ]
lp_imaging<-lp_imaging[!grepl( "_Manders_" , rownames(lp_imaging)) ,  ]
lp_imaging<-lp_imaging[!grepl( "_RWC_" , rownames(lp_imaging)) ,  ]
lp_imaging<-lp_imaging[!grepl( "SmallBODIPY" , rownames(lp_imaging) ) ,  ]

mean<-rowMeans(lp_imaging)
var<-cbind(mean, lp_imaging)
lp_imaging<-subset(var, var$mean != "0")

lp_imaging<- lp_imaging[complete.cases(lp_imaging), ]
writeLines(sprintf("lp_imaging; [%d, %d ]", dim(lp_imaging)[1], dim(lp_imaging)[2]))
lp_imaging<-as.data.frame(t(lp_imaging[,2:dim(lp_imaging)[2]]))
curr_lp_data<-cbind(lp_meta, lp_imaging)
curr_lp_data <- curr_lp_data[complete.cases(curr_lp_data), ] 

### iterate on conditions ###
days <- c(0, 3, 8, 14)
depots <- c("vc", "sc")
sex_stratifications <- c(TRUE, FALSE) 
FFA_tags <- c(0, 1)

day <- 0
depot <- "sc"
FFA_tag <- 0
for (day in days){
  for (depot in depots){
    for (sex_stratified in sex_stratifications){
      for (FFA_tag in FFA_tags){
        current_model_input<-subset(curr_lp_data, curr_lp_data$cellType == depot & curr_lp_data$Day == day & curr_lp_data$FFA == FFA_tag)
        if (sex_stratified){
          current_model_input_male <- subset(current_model_input, current_model_input$sex == 1)
          current_model_input_female <- subset(current_model_input, current_model_input$sex == 2)
          if (dim(current_model_input_male)[1] > 5){
            writeLines(sprintf("D%s-%s-Male-%s (N=%d)", day, depot, FFA_tag, dim(current_model_input_male)[1]))
            run_regression(current_model_input_male, gene_feature_output_path, depot, day, "Male", FFA_tag, gene_ensb_id)            
          }
          if (dim(current_model_input_female)[1] > 5){
            writeLines(sprintf("D%s-%s-Female-%s (N=%d)", day, depot, FFA_tag, dim(current_model_input_female)[1]))
            run_regression(current_model_input_female, gene_feature_output_path, depot, day, "Female", FFA_tag, gene_ensb_id)
          }
        }else{
          # no sex stratification, use the current_model_input
          if (dim(current_model_input)[1] > 5){
            writeLines(sprintf("D%s-%s-both-%s (N=%d)", day, depot, FFA_tag, dim(current_model_input)[1]))
            run_regression(current_model_input, gene_feature_output_path, depot, day, "both_sex", FFA_tag, gene_ensb_id)            
          }
        }        
      }
    }
  }
}  

##### PLOT GENE-FEATURE RESULTS ####
### iterate on conditions ###
days <- c(0, 3, 8, 14)
depots <- c("vc", "sc")
sex_stratifications <- c(TRUE, FALSE) 
FFA_tags <- c(0, 1)

for (day in days){
  for (depot in depots){
    for (FFA_tag in FFA_tags){
      for (sex_label in c("Male", "Female", "both_sex")){
        file_path <- sprintf("%s/%s_%d_%s_%s.csv", gene_feature_output_path, depot, day, sex_label, FFA_tag)
        if(file.exists(file_path)){
          results_df <- read.csv(file_path)
          significant_results_df <- results_df[which(results_df$q_value < .05),]
          
          if (dim(significant_results_df)[1] > 5){
            out_file_path <- sprintf("%s/%s_D%d_%s_%s", gene_feature_output_path, depot, day, sex_label, FFA_tag)
            in_title <- sprintf("%s_D%d_%s_%s", depot, day, sex_label, FFA_tag)
            MyPieDonut(significant_results_df, in_title, out_file_path)
          }
          pos <- significant_results_df[significant_results_df$beta>0,]
          if (dim(pos)[1] > 5){
            out_file_path <- sprintf("%s/%s_D%d_%s_%s_pos", gene_feature_output_path, depot, day, sex_label, FFA_tag)
            in_title <- sprintf("%s_D%d_%s_%s", depot, day, sex_label, FFA_tag)
            MyPieDonut(pos, in_title, out_file_path)
          }
          neg <- significant_results_df[significant_results_df$beta<=0,]
          if (dim(neg)[1] > 5){
            out_file_path <- sprintf("%s/%s_D%d_%s_%s_neg", gene_feature_output_path, depot, day, sex_label, FFA_tag)
            in_title <- sprintf("%s_D%d_%s_%s", depot, day, sex_label, FFA_tag)
            MyPieDonut(neg, in_title, out_file_path)
            write.csv(neg, sprintf("%s_resultdf.csv", out_file_path))
          }
        }
      }
    }
  }
}
