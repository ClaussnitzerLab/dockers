library(sjstats)
library(cowplot)
library(dplyr)
library(ggplot2)
library(rlang)
library(qvalue)
library(ggsignif)
library(stringr)

rm(list = ls())

##### volcano plaots ######
plot_volcano <- function(volcano, image_path, image_title){
  ####define color categories###
  volcano <- volcano %>% dplyr::mutate(feature_color = 
                                         ifelse(grepl("BODIPY", variable),'Lipid', 
                                                ifelse(grepl("AGP", variable),'AGP', 
                                                       ifelse(grepl("Mito", variable),'Mito',
                                                              ifelse(grepl("DNA", variable),'DNA',
                                                                     'other')))))
  volcano$feature_color<-factor(volcano$feature_color, levels = c("Mito","AGP","Lipid","DNA", "other"))
  ####remove missing values###
  volcano<- volcano[complete.cases(volcano), ]
  ####set significance level###
  
  volcano <- volcano %>% dplyr::mutate(feature_sig = ifelse(volcano$q_value <= 0.05,'5%FDR','not'))
  ll <- subset(volcano, volcano$feature_sig == '5%FDR')
  plot <- ggplot(volcano, aes(x=beta, y=-log10(q_value))) +  # `t-test` 
    geom_point(aes(color=feature_color,
                   size = feature_sig,
                   alpha = feature_sig))+
    
    xlab("estimate") +
    ylab("-log10 Adj p-value") +
    #ylim(0,4)+
    #xlim(-5,5)+
    scale_size_manual(name = "", 
                      values = c("5%FDR" = 1.5, "not" = 1))+
    scale_alpha_manual(name = "", 
                       values = c("5%FDR" = 0.75, "not" = 0.1))+
    scale_color_manual(name = "", 
                       values = c("Mito" = "#f56464" ,
                                  "AGP" = "#f2cf41",
                                  "Lipid" = "#4dac26",
                                  "DNA" = "#65a6db",
                                  "other" = "#959ca3"))+
    theme_bw() +
    theme(strip.background = element_rect(colour = "black", fill = "#fdfff4"))
  # sprintf("Depot: %s, D%d, sex: %s, %s vs. %s, N=%d, %s", depot, day, sex_option, group_class1, group_class2, dim(data_anov)[1], header)
  plot+ggtitle(image_title)
  ggsave(image_path)  
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

# results_df <- significant_results_df, title, 
# file_name <- out_file_path
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
####### End of DONUT PLOT #######

#### generating images for a set of identified features
generate_results <- function(features, output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header){
  complete_features_path <- sprintf("%s/tables/features_%s_%s_%s_FFA%s_%s_%s.csv", output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple)
  write.csv(features, complete_features_path)
  
  image_path <- sprintf("%s/figures/features_%s_%s_%s_FFA%s_%s_%s_%s.pdf", output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
  image_title <- sprintf("Depot:%s;D%s;%s;FFA%s;%s vs. %s;%s", depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
  plot_volcano(features, image_path, image_title)
  
  significant_features <- features[which(features$q_value <= .05),]
  if (dim(significant_features)[1] > 5){
    out_file_path <- sprintf("%s/figures/features_%s_%s_%s_FFA%s_%s_%s_donut_%s", output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    in_title <- sprintf("%sD%s;%s;FFA%s;%s_%s_donut_%s", depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    MyPieDonut(significant_features, in_title, out_file_path)
  }
  pos <- significant_features[significant_features$beta>0,]
  if (dim(pos)[1] > 5){
    out_file_path <- sprintf("%s/figures/features_%s_%s_%s_FFA%s_%s_%s_donut_positive_%s", output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    in_title <- sprintf("%sD%s;%s;FFA%s;Pos;%s_%s_donut_%s", depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    MyPieDonut(pos, in_title, out_file_path)
    write.csv(pos, sprintf("%s_resultdf.csv", out_file_path))
  }
  neg <- significant_features[significant_features$beta<=0,]
  if (dim(neg)[1] > 5){
    out_file_path <- sprintf("%s/figures/features_%s_%s_%s_FFA%s_%s_%s_donut_negative_%s", output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    in_title <- sprintf("%sD%s;%s;FFA%s;Neg;%s_%s_donut_%s", depot, day, sex_option, FFA_tag, g1_simple, g2_simple, header)
    MyPieDonut(neg, in_title, out_file_path)
    write.csv(neg, sprintf("%s_resultdf.csv", out_file_path))
  }  
}


update_regression_arrays <- function(current_model_input, outcome, sex_label){
  current_model_input$sex <- as.factor(current_model_input$sex)
  # current_model_input$T2D <- as.factor(current_model_input$T2D)
  current_model_input$batch <- as.factor(current_model_input$batch)
  current_model_input <- na.omit(current_model_input)
  formula <- sprintf("get('%s') ~ genotype + age + BMI", outcome)
  if (sex_label == "both"){  #Changed here
    if (length(unique(current_model_input$sex)) > 1){
      formula <- sprintf("%s + sex", formula)
    }        
  }
  if (length(unique(current_model_input$batch)) > 1){
    formula <- sprintf("%s + batch", formula)
  } 
  lm_model <- summary.lm(lm(formula = formula, data = current_model_input))
  beta <-lm_model$coefficients[2, 1]
  pval <-lm_model$coefficients[2, 4]
  
  output <- list(beta, 0, pval, FALSE)
  return(output)
}

# runnning lm 
# curr_data <- combined_allelees, 
# group_class1 <- group01, 
# group_class2 <- group11, 
# sex_label <- sex_option
run_regression <- function(curr_data, group_class1, group_class2, sex_label){
  case <- subset(curr_data, curr_data$genotype == group_class1)
  control <- subset(curr_data, curr_data$genotype == group_class2)
  output_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(output_df) <- c("variable", "coeff", "pvalue")
  if (dim(case)[1] > 1 & dim(control)[1] > 1){
    data_anov<-rbind(case, control)
    #data_anov$genotype<-factor(data_anov$genotype, levels= c(group_class1, group_class2))
    # imaging_data <- select(data_anov, -c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression", "genotype"))
    imaging_data <- select(data_anov, -c("batch", "newcol", "patientID", "sex", "age", "BMI", "FFA", "cellType", "Day", "expression", "genotype"))
    for (outcome in colnames(imaging_data)){
      output <- update_regression_arrays(data_anov, outcome, sex_label)
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
    p<-output_df$pval
    if (is_empty(p)){
      model_output_final <- output_df
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
    output_df <- model_output_final
  }
  return(output_df)
}

# run analysis on a set of features
# curr_data, output_path, depot, day, sex_option, FFA_tag
run_on_a_set <- function(curr_data, output_path, depot, day, sex_option, FFA_tag){
  row <- 1
  for (row in 1:nrow(not_combined_group_definitions)) {
    g1 <- not_combined_group_definitions[row, "Group0"]
    g2  <- not_combined_group_definitions[row, "Group1"]
    g1_simple <- not_combined_group_definitions[row, "SimpGroup0"]
    g2_simple <- not_combined_group_definitions[row, "SimpGroup1"]
    features <- run_regression(curr_data, g1, g2, sex_option)
    if (dim(features)[1] > 0){
      generate_results(features, output_path, depot, day, sex_option, FFA_tag, g1_simple, g2_simple, "")
    }
  }
  
  combined_allelees <- curr_data
  combined_allelees["genotype"] <- ifelse(combined_allelees$genotype == group00, group01, combined_allelees$genotype)
  header <- sprintf("combined_%s_%s", simple_group00, simple_group01)
  features <- run_regression(combined_allelees, group01, group11, sex_option)
  if(dim(features)[1] > 0){
    generate_results(features, output_path, depot, day, sex_option, FFA_tag, sprintf("combined_%s_%s", simple_group00, simple_group01), simple_group11, header)
  }
  
  
  combined_allelees <- curr_data
  combined_allelees["genotype"] <- ifelse(combined_allelees$genotype == group11, group01, combined_allelees$genotype)
  header <- sprintf("combined_%s_%s", simple_group11, simple_group01)
  features <- run_regression(combined_allelees, group01, group00, sex_option)
  if(dim(features)[1] > 0){
    generate_results(features, output_path, depot, day, sex_option, FFA_tag, sprintf("combined_%s_%s", simple_group11, simple_group01), simple_group00, header)
  }
}

# 
run_analsys <- function(gene_name, gene_ensb_id, snp_name, allele, output_path, snp_coordinate, snp_codes_00, snp_codes_11){
  allele["ID"] <- rownames(allele)
  l <- merge(data, allele, by.x = "patientID", by.y = "ID")
  l["expression"] <- l[gene_ensb_id]
  summary(l["expression"])
  
  l <- l %>% select(-contains("ENSG"))
  
  ll0 <- snp_coordinate
  ll1 <- as.data.frame(summary(as.factor(l$genotype)))
  l["genotype"] <- ifelse(l$genotype=="0|0", sprintf("%s/%s", snp_codes_00, snp_codes_00), 
                          ifelse(l$genotype=="1|1", sprintf("%s/%s", snp_codes_11, snp_codes_11), 
                                 sprintf("%s/%s", snp_codes_00, snp_codes_11)))
  ll2 <- as.data.frame(summary(as.factor(l$genotype)))
  summary(as.factor(l$genotype))
  
  #lp_meta <- select(l, c("batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression", "genotype"))
  #lp_imaging <- select(l, -c("X","batch", "newcol", "patientID", "sex", "age", "BMI", "T2D", "FFA", "cellType", "Day", "expression", "genotype"))
  lp_meta <- select(l, c("batch", "newcol", "patientID", "sex", "age", "BMI", "FFA", "cellType", "Day", "expression", "genotype"))
  lp_imaging <- select(l, -c("X","batch", "newcol", "patientID", "sex", "age", "BMI", "FFA", "cellType", "Day", "expression", "genotype"))
  
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
  
  lp_imaging<-as.data.frame(t(lp_imaging[,2:dim(lp_imaging)[2]]))
  curr_lp_data<-cbind(lp_meta, lp_imaging)
  curr_lp_data <- curr_lp_data[complete.cases(curr_lp_data), ]
  
  curr_lp_data$sex <- ifelse(curr_lp_data$sex == 2, "Female", "Male")
  output_table <- as.data.frame(t(rbind(ll0, ll1, ll2)))
  
  write.csv(output_table, sprintf("%s/allele_info_imaging.csv", output_path), row.names = F) 
  write.csv(curr_lp_data, sprintf("%s/genotyping_data_imaging.csv", output_path))  
  
  days <- c(0, 3, 8, 14)
  depots <- c("vc", "sc")
  sex_stratifications <- c(TRUE, FALSE) 
  FFA_tags <- c(0, 1)
  
  
  
  day <- 14
  depot <- "sc"
  FFA_tag <- 0
  sex_stratified <- FALSE
  for (day in days){
    for (depot in depots){
      for (sex_stratified in sex_stratifications){
        for (FFA_tag in FFA_tags){
          curr_data <- subset(curr_lp_data, curr_lp_data$Day == day & curr_lp_data$cellType == depot & curr_lp_data$FFA == FFA_tag)
          if (!sex_stratified){  # both sex
            sex_option <- "both"
            run_on_a_set(curr_data, output_path, depot, day, sex_option, FFA_tag)
            
          }else{  # sex stratifies
            sex_option <- "Male"
            curr_data_male <- subset(curr_data, curr_data$sex == "Male")
            run_on_a_set(curr_data_male, output_path, depot, day, sex_option, FFA_tag)
            sex_option <- "Female"
            curr_data_Female <- subset(curr_data, curr_data$sex == "Female")
            run_on_a_set(curr_data_Female, output_path, depot, day, sex_option, FFA_tag)
          }
        }
      }
    }
  }
}

#curr_plot_data <- curr_temp
# combined_allelees, snp_name, image_title, image_path, day, depot, FFA_tag, expression_result_table
plot_expression_per_allelee <- function(curr_plot_data, snp_name, image_title, image_path, day, depot, FFA_tag, expression_result_table){
  curr_plot_data["genotype"] <- as.factor(curr_plot_data$genotype)
  # summary(as.factor(curr_plot_data$genotype))
  # writeLines(sprintf("%d", length(unique(curr_plot_data$genotype))))
  base <- max(curr_plot_data$expression)
  
  if (length(unique(curr_plot_data$genotype)) == 2){
    res <- wilcox.test(expression ~ genotype, data = curr_plot_data)
    temp <- data.frame(day, depot, FFA_tag, image_title, res$p.value)
    colnames(temp) <- colnames(expression_result_table)
    expression_result_table <- rbind(expression_result_table, temp)
  }
  ggplot(curr_plot_data, aes(x=genotype, y=expression)) + 
    geom_boxplot() + ggtitle(image_title) + 
    xlab("Genotype") + ylab(sprintf("Expression %s(%s;%s)", gene_name, gene_ensb_id, snp_name)) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
    geom_signif(comparisons=list(c(sprintf("%s/%s", snp_codes_00, snp_codes_00), sprintf("%s/%s", snp_codes_11, snp_codes_11))),  tip_length = 0, vjust=0.1, y_position = base) + 
    geom_signif(comparisons=list(c(sprintf("%s/%s", snp_codes_00, snp_codes_00), sprintf("%s/%s", snp_codes_00, snp_codes_11))),  tip_length = 0, vjust=0.1, y_position = base+floor(base/3)) +
    geom_signif(comparisons=list(c(sprintf("%s/%s", snp_codes_11, snp_codes_11), sprintf("%s/%s", snp_codes_00, snp_codes_11))),  tip_length = 0, vjust=0.1, y_position = base+floor(base/2))
  ggsave(image_path)
  return(expression_result_table)
}

clean_db <- function(temp){
  l <- unique(as.factor(temp$genotype))
  to_be_removed <- c()
  for (LL in l){
    if (dim(temp[temp$genotype == LL,])[1] < 2){
      to_be_removed <- c(to_be_removed, LL)
    }
  }
  for (c in to_be_removed){
    temp <- temp[temp$genotype != c,]
  }
  return(temp)
}

plot_eqtl <- function(curr_data, output_path, depot, day, sex_option, FFA_tag, snp_name, snp_codes_00, snp_codes_11, expression_result_table, formula){
  #curr_data <- curr_data_male, output_path, depot, day, sex_option, FFA_tag, snp_name, snp_codes_00, snp_codes_11, expression_result_table
  if (length(unique(curr_data$sex))>1){
    formula <- sprintf("%s + sex", formula)
  }
  model <- lm(eval(formula), data = curr_data)
  curr_data$expression <- curr_data$expression-predict(model, newdata = curr_data)
  
  curr_temp <- clean_db(curr_data)
  if (length(unique(curr_temp$genotype )) >= 2){
    image_title <- sprintf("%s;%s;%s;%s", depot, day, sex_option, FFA_tag)
    image_path <- sprintf("%s/figures/%s_%s_%s_%s.pdf", output_path, depot, day, sex_option, FFA_tag)
    expression_result_table <- plot_expression_per_allelee(curr_data, snp_name, image_title, image_path, day, depot, FFA_tag, expression_result_table)
  }
  
  #### combining minor allelee and heter
  major <- sprintf("%s/%s", snp_codes_00, snp_codes_00)
  minor <- sprintf("%s/%s", snp_codes_11, snp_codes_11)
  heter <- sprintf("%s/%s", snp_codes_00, snp_codes_11)
  
  combined_allelees <- curr_data
  combined_allelees["genotype"] <- ifelse(combined_allelees$genotype == minor, heter, combined_allelees$genotype)
  header <- sprintf("%s merged with %s", minor, heter)
  if (length(unique(combined_allelees$genotype )) >= 2){
    image_title <- sprintf("%s;%s;%s;%s;%s", depot, day, sex_option, FFA_tag, header)
    
    image_path <- sprintf("%s/figures/%s_%s_%s_%s_%s.pdf", output_path, depot, day, sex_option, FFA_tag, sprintf("combined_%s_%s", 
                                                                                                                 str_replace(minor, "/", ""), 
                                                                                                                 str_replace(heter, "/", "")))
    expression_result_table <- plot_expression_per_allelee(combined_allelees, snp_name, image_title, image_path, day, depot, FFA_tag, expression_result_table)
  }
  
  combined_allelees <- curr_data
  combined_allelees["genotype"] <- ifelse(combined_allelees$genotype == major, heter, combined_allelees$genotype)
  header <- sprintf("%s merged with %s", major, heter)
  if (length(unique(combined_allelees$genotype )) >= 2){
    image_title <- sprintf("%s;%s;%s;%s;%s", depot, day, sex_option, FFA_tag, header)
    image_path <- sprintf("%s/figures/%s_%s_%s_%s_%s.pdf", output_path, depot, day, sex_option, FFA_tag, sprintf("combined_%s_%s", 
                                                                                                                 str_replace(major, "/", ""), 
                                                                                                                 str_replace(heter, "/", "")))
    expression_result_table <- plot_expression_per_allelee(combined_allelees, snp_name, image_title, image_path, day, depot, FFA_tag, expression_result_table)
  }
  return(expression_result_table)
}

run_eqtl_analysis <- function(gene_name, gene_ensb_id, snp_name, allele, output_path, snp_coordinate, snp_codes_00, snp_codes_11, expression_result_table){
  allele["ID"] <- rownames(allele)
  # names(allele)[names(allele) == snp_coordinate] <- "genotype"
  l <- merge(data, allele, by.x = "patientID", by.y = "ID")
  l["expression"] <- l[gene_ensb_id]
  
  ll0 <- snp_coordinate
  ll1 <- as.data.frame(summary(as.factor(l$genotype)))
  l["genotype"] <- ifelse(l$genotype=="0|0", sprintf("%s/%s", snp_codes_00, snp_codes_00), 
                          ifelse(l$genotype=="1|1", sprintf("%s/%s", snp_codes_11, snp_codes_11), 
                                 sprintf("%s/%s", snp_codes_00, snp_codes_11)))
  ll2 <- as.data.frame(summary(as.factor(l$genotype)))
  summary(as.factor(l$genotype))
  
  curr_genotype_data <- select(l, c("genotype", "expression", "FFA", "cellType", "Day", "sex", "patientID", "age", "BMI", "batch"))
  curr_genotype_data$sex <- ifelse(curr_genotype_data$sex == 2, "Female", "Male")
  output_table <- as.data.frame(t(rbind(ll0, ll1, ll2)))
  
  write.csv(output_table, sprintf("%s/allele_info.csv", output_path))
  write.csv(curr_genotype_data, sprintf("%s/genotyping_data.csv", output_path))
  
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
          if (FFA_tag == 1 & day == 0) next;
          curr_data <- subset(curr_genotype_data, curr_genotype_data$Day == day & curr_genotype_data$cellType == depot & curr_genotype_data$FFA == FFA_tag)
          # this section adjusts for age, bmi, batch
          curr_data$age <- as.numeric(curr_data$age)
          curr_data$BMI <- as.numeric(curr_data$BMI)
          curr_data$batch <- as.factor(curr_data$batch)
          curr_data <- curr_data[complete.cases(curr_data),]
          
          formula <- "expression ~ age + BMI"
          if (length(unique(curr_data$batch)) > 1){
            formula <- sprintf("%s + batch", formula)
          }

          if (!sex_stratified){  # both sex
            sex_option <- "both"
            expression_result_table <- plot_eqtl(curr_data, output_path, depot, day, sex_option, FFA_tag, snp_name, snp_codes_00, snp_codes_11, expression_result_table, formula)
          }else{  # sex stratifies
            sex_option <- "Male"
            curr_data_male <- subset(curr_data, curr_data$sex == "Male")
            expression_result_table <-  plot_eqtl(curr_data_male, output_path, depot, day, sex_option, FFA_tag, snp_name, snp_codes_00, snp_codes_11, expression_result_table, formula)
            sex_option <- "Female"
            curr_data_Female <- subset(curr_data, curr_data$sex == "Female")
            expression_result_table <- plot_eqtl(curr_data_Female, output_path, depot, day, sex_option, FFA_tag, snp_name, snp_codes_00, snp_codes_11, expression_result_table, formula)
          }
        }
      }
    }
  }
  return(expression_result_table)
}


args = commandArgs(trailingOnly=TRUE)
LP_data_path <- args[1]
input_vcf_file <- args[2]
gene_name <- args[3]
snp_name <-  args[4]
gene_ensb_id <- args[5]
snp_coordinate <- args[6]
snp_codes_00 <- args[7]
snp_codes_11 <- args[8]



root_output <- "LP_output_local"
dir.create(root_output)

data<-read.csv(LP_data_path)
summary(as.factor(data$sex))
output_path <- sprintf("%s/snp_features_%s", root_output, snp_name)
dir.create(file.path(output_path))

group00 <- sprintf("%s/%s", snp_codes_00, snp_codes_00)
group11 <- sprintf("%s/%s", snp_codes_11, snp_codes_11)
group01 <- sprintf("%s/%s", snp_codes_00, snp_codes_11)

simple_group00 <- sprintf("%s%s", snp_codes_00, snp_codes_00)
simple_group11 <- sprintf("%s%s", snp_codes_11, snp_codes_11)
simple_group01 <- sprintf("%s%s", snp_codes_00, snp_codes_11)

Group0 <- c(group00, group00, group11)
Group1 <- c(group01, group11, group01)
SimpGroup0 <- c(simple_group00, simple_group00, simple_group11)
SimpGroup1 <- c(simple_group01, simple_group11, simple_group01)
not_combined_group_definitions <- data.frame(Group0, Group1, SimpGroup0, SimpGroup1)

df <- select(data, c("sex", "Day", "cellType", gene_ensb_id))
colnames(df) <- c("sex", "Day", "cellType", "expression")
df$Day <- as.factor(df$Day)
df$sex <- ifelse(df$sex == 2, "Female", "Male")
df$sex <- as.factor(df$sex)
ggplot(df, aes(x=Day, y=expression, fill=sex)) + geom_boxplot()+facet_grid(~cellType)
ggsave(sprintf("%s/expression_%s.pdf", output_path, gene_ensb_id))


vcf <- read.csv(input_vcf_file)
rownames(vcf) <- vcf$IDs
allele <- vcf # as.data.frame(vcfR2loci(vcf))
dir.create(file.path(sprintf("%s/figures", output_path)))
dir.create(file.path(sprintf("%s/tables", output_path)))
run_analsys(gene_name, gene_ensb_id, snp_name, allele, output_path, snp_coordinate, snp_codes_00, snp_codes_11)
run_eqtl_analysis(gene_name, gene_ensb_id, snp_name, allele, output_path, snp_coordinate,   snp_codes_00, snp_codes_11, expression_result_table)
write.csv(expression, sprintf("%s/eqtl_p.csv", output_path), row.names = F)






write.csv(expression, sprintf("%s/eqtl_p.csv", output_path), row.names = F)


