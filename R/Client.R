#' @title DisTxRESP Client
#' @description An R6 client to interact with the DisTxRESP database. It allows users to query metadata, filter datasets, batch-download files, and load data directly into R memory for downstream visualization and analysis.
#' @importFrom R6 R6Class
#' @importFrom data.table fread
#' @importFrom base64enc base64decode
#' @importFrom httr GET http_error http_status progress
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw labs geom_vline geom_hline geom_boxplot geom_jitter theme_classic
#' @importFrom ggpubr stat_compare_means
#' @importFrom pROC roc ggroc
#' @importFrom rlang .data
#' @export
DisTxRESP <- R6Class("DisTxRESP",
                     private = list(
                       api_token = "aHR0cDovLzIxOC44LjI0MS4yNDg6MzgzOC9EaXNUeFJFU1AvZG93bmxvYWQv",
                       make_request = function(endpoint) {
                         dl_url = paste0(rawToChar(base64decode(private$api_token)), endpoint)
                         return(dl_url)
                       },
                       download_file = function(url, destfile) {
                         response <- httr::GET(url,
                                               httr::write_disk(destfile, overwrite = TRUE),
                                               httr::progress(),
                                               httr::config(max_recv_speed_large = 1000000))
                         if(httr::http_error(response)) {
                           stop("Download failed: ", httr::http_status(response)$message)
                         }
                         return(invisible(response))
                       }
                     ),
                     public = list(
                       #' @field metadata A data.frame containing the full metadata manifest downloaded from the server.
                       metadata = NULL,
                       #' @field sub_table A data.frame containing the filtered subset of the metadata.
                       sub_table = NULL,
                       #' @field local_data A nested list storing the loaded clinical data, feature matrices, and differential results in R memory.
                       local_data = NULL,

                       #' @description Create a new database client and the fetch full metadata.
                       #' @param manifest String. The local file path for the metadata manifest. If the file already exists locally, it will be read directly; otherwise, it will be downloaded from the server. Defaults to "Download_Manifest.csv".
                       #' @return A new \code{DisTxRESP} object.
                       initialize = function(manifest="Download_Manifest.csv") {
                         if (!file.exists(manifest)) {
                           message("Connecting to server and fetching metadata table...")
                           meta_url <- private$make_request("Download_Manifest.csv")

                           download_file(url = meta_url, destfile = manifest)

                         }else{
                           message("Local file found, reading...")
                         }
                         self$metadata <- fread(manifest,data.table = FALSE,header = TRUE)
                         self$sub_table <- self$metadata
                       },

                       #' @description Filter the metadata table based on specific biological or clinical requirements.
                       #' @param dtxr Character vector. The DTXR dataset identifiers from the DisTxRESP database (e.g., c("DTXR100001")).
                       #' @param data_type Character vector. The omics data types to retain (e.g., c("Transcriptomics")"").
                       #' @param disease Character vector. The primary disease types to filter by (e.g., c("Leukemia")).
                       #' @param disease_sub Character vector. The specific disease subtypes (e.g., c("Acute myeloid leukemia")).
                       #' @param treatment Character vector. The treatment regimens applied to the cohorts (e.g., c("Acute myeloid leukemia")).
                       #' @param intervention Character vector. The intervention types (e.g., c("Targeted Therapy")).
                       #' @param sampling_location Character vector. The tissue origins or sampling locations(e.g., c("Tissue","Bone Marrow")).
                       #' @param feature Character vector. The feature types to retain (e.g., c("Gene Expression","Clinical Data")).
                       #' @param file_type Character vector. The file types to retain (e.g., c("Meta Info", "Feature Matrix")).
                       #' @param min_size Integer. The minimum sample size required for a dataset to be retained. Defaults to 0.
                       #' @return Returns the modified R6 object invisibly, allowing for method chaining.
                       filter_metadata = function(dtxr=NULL,data_type = NULL,disease = NULL,
                                                  disease_sub = NULL,treatment = NULL,
                                                  intervention = NULL,sampling_location = NULL,
                                                  feature = NULL, file_type = NULL,min_size=0) {
                         if (is.null(self$metadata)) stop("Metadata is empty. Please re-initialize.")
                         temp_df <- self$metadata

                         if (!is.null(dtxr)) {
                           temp_df <- temp_df[temp_df$Dataset %in% dtxr, ]
                         }
                         if (!is.null(data_type)) {
                           temp_df <- temp_df[grepl(paste0(data_type,collapse = "|"),temp_df$Omics), ]
                         }
                         if (!is.null(disease)) {
                           temp_df <- temp_df[grepl(paste0(disease,collapse = "|"),temp_df$Disease_Type), ]
                         }
                         if (!is.null(disease_sub)) {
                           temp_df <- temp_df[grepl(paste0(disease_sub,collapse = "|"),temp_df$Disease_Subtype), ]
                         }
                         if (!is.null(treatment)) {
                           temp_df <- temp_df[grepl(paste0(treatment,collapse = "|"),temp_df$Treatment_Regimen), ]
                         }
                         if (!is.null(intervention)) {
                           temp_df <- temp_df[grepl(paste0(intervention,collapse = "|"),temp_df$Intervention), ]
                         }
                         if (!is.null(sampling_location)) {
                           temp_df <- temp_df[grepl(paste0(sampling_location,collapse = "|"),temp_df$Sampling_Location), ]
                         }
                         if (!is.null(feature)) {
                           temp_df <- temp_df[grepl(paste0(feature,collapse = "|"),temp_df$Feature_Type), ]
                         }
                         if (!is.null(file_type)) {
                           temp_df <- temp_df[grepl(paste0(file_type,collapse = "|"),temp_df$File_Type), ]
                         }
                         temp_df = temp_df[temp_df$Sample_Size>=min_size,]
                         # if (!is.null(dtxr)) {
                         #   private$check_fuzzy_match(dtxr, temp_df$Dataset, "Dataset (dtxr)")
                         #   temp_df <- temp_df[temp_df$Dataset %in% dtxr, ]
                         # }
                         # if (!is.null(data_type)) {
                         #   private$check_fuzzy_match(data_type, temp_df$Omics, "Omics (data_type)")
                         #   temp_df <- temp_df[temp_df$Omics %in% data_type, ]
                         # }
                         # if (!is.null(disease)) {
                         #   private$check_fuzzy_match(disease, temp_df$Disease_Type, "Disease_Type (disease)")
                         #   temp_df <- temp_df[temp_df$Disease_Type %in% disease, ]
                         # }
                         # if (!is.null(feature)) {
                         #   private$check_fuzzy_match(feature, temp_df$Feature_Type, "Feature_Type (feature)")
                         #   temp_df <- temp_df[temp_df$Feature_Type %in% feature, ]
                         # }
                         # if (!is.null(file_type)) {
                         #   private$check_fuzzy_match(file_type, temp_df$File_Type, "File_Type (file_type)")
                         #   temp_df <- temp_df[temp_df$File_Type %in% file_type, ]
                         # }

                         self$sub_table <- temp_df
                         message(sprintf("Filtered down to %d records.", nrow(self$sub_table)))

                         return(invisible(self))
                       },

                       #' @description Batch download files in the filtered \code{sub_table} to a local directory.
                       #' @param output_dir String. The target local directory path to save the downloaded files. Defaults to "downloaded_data".
                       #' @param skip_existing Logical. Whether to skip downloading files that already exist in the target directory. Defaults to TRUE.
                       #' @return Returns the modified R6 object invisibly.
                       download_files = function(output_dir = "downloaded_data",skip_existing=TRUE) {
                         if (nrow(self$sub_table) == 0) stop("The filtered table is empty. Nothing to download.")

                         for (i in seq_len(nrow(self$sub_table))) {
                           row_data <- self$sub_table[i, ]
                           data_id <- row_data$Dataset
                           f_type <- gsub(" ","_",row_data$File_Type)
                           f_name <- row_data$File_Name

                           temp_dir <- ifelse(f_type=="Meta_Info",paste0(output_dir,"/",data_id),paste0(output_dir,"/",data_id,"/",f_type))
                           file_name <- paste0(temp_dir,"/",f_name)
                           if (skip_existing) {
                             if (file.exists(file_name)) {
                               message(sprintf("File already exists: %s...", f_name))
                               next
                             }
                           }

                           if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

                           dl_url <- private$make_request(sprintf("%s/%s",data_id, f_name))

                           message(sprintf("Downloading %s...", f_name))
                           private$download_file(url = dl_url,destfile = file_name)
                         }
                         message(" \n Batch download complete.")
                         return(invisible(self))
                       },

                       #' @description Fetch the filtered files and load them into an R list structure.
                       #' @param local_dir String. The directory path containing locally downloaded files.
                       #' @return Returns the modified R6 object invisibly. The loaded data is accessible via the \code{local_data} field.
                       load_to_memory = function(local_dir="downloaded_data") {
                         if (nrow(self$sub_table) == 0) stop("The filtered table is empty. Nothing to load.")

                         result_list <- list()

                         for (i in seq_len(nrow(self$sub_table))) {
                           row_data <- self$sub_table[i, ]
                           data_id <- row_data$Dataset
                           f_type <- gsub(" ","_",row_data$File_Type)
                           f_name <- row_data$File_Name
                           feature <- row_data$Feature_Type

                           if (!is.null(local_dir)) {
                             temp_dir <- ifelse(f_type=="Meta_Info",paste0(local_dir,"/",data_id),paste0(local_dir,"/",data_id,"/",f_type))
                             file_name <- paste0(temp_dir,"/",f_name)
                             if (file.exists(file_name)) {
                               message(sprintf("File found, loading %s into memory...", f_name))
                             }else{
                               stop(sprintf("The file was not found: %s", f_name))
                             }
                           }else{
                             message("Loading path cannot be empty!")
                           }

                           tryCatch({
                             if (f_type=="Meta_Info") {
                               result_list[[data_id]][[f_type]] <- fread(file_name, data.table = F, header = T)
                             }else{
                               result_list[[data_id]][[f_type]][[feature]] <- fread(file_name, data.table = F, header = T)
                             }
                           }, error = function(e) {
                             warning(sprintf("Failed to load %s: %s", data_id, e$message))
                           })
                         }
                         message(" \n All selected files successfully loaded into memory.")
                         self$local_data = result_list
                         return(invisible(self))
                       },

                       #' @description Generate a Volcano Plot for differential analysis results.
                       #' @param dataset String. The dataset identifier (e.g., "DTXR100001").
                       #' @param feature_type String. The type of feature (e.g., "Gene Expression").
                       #' @param logX_col String. The column name representing the log2 Fold Change (or LogOR) in the differential matrix. Defaults to "LogOR".
                       #' @param pval_col String. The column name representing the P-value or FDR in the differential matrix. Defaults to "P.value".
                       #' @param fc_cutoff Numeric. The absolute threshold for fold change significance. Defaults to 1.0.
                       #' @param p_cutoff Numeric. The threshold for P-value significance. Defaults to 0.05.
                       #' @return A \code{ggplot} object.
                       plot_volcano = function(dataset, feature_type, logX_col = "LogOR", pval_col = "P.value",
                                               fc_cutoff = 1.0, p_cutoff = 0.05) {

                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         diff_df <- self$local_data[[dataset]]$Differential_Results[[feature_type]]
                         if (is.null(diff_df)) stop("No difference matrix found for this feature type!")

                         diff_df$Significance <- "Not Significant"
                         diff_df$Significance[diff_df[[logX_col]] > fc_cutoff & diff_df[[pval_col]] < p_cutoff] <- "Up"
                         diff_df$Significance[diff_df[[logX_col]] < -fc_cutoff & diff_df[[pval_col]] < p_cutoff] <- "Down"

                         p <- ggplot2::ggplot(diff_df, ggplot2::aes(x = .data[[logX_col]],
                                                                    y = -log10(.data[[pval_col]]),
                                                                    color = .data[["Significance"]])) +
                           ggplot2::geom_point(alpha = 0.8, size = 1.5) +
                           ggplot2::scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "Not Significant" = "grey80")) +
                           ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
                           ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
                           ggplot2::theme_bw() +
                           ggplot2::labs(title = sprintf("%s (%s)", dataset, feature_type),
                                         x = logX_col, y = "-log10(P-value)")

                         return(p)
                       },

                       #' @description Generate a Boxplot to visualize feature expression across different clinical groups.
                       #' @param dataset String. The dataset identifier (e.g., "DTXR100001").
                       #' @param feature_type String. The type of feature matrix (e.g., "Gene Expression").
                       #' @param feature_name String. The specific feature name to plot (e.g., "CD274").
                       #' @param group_col String. The column name in the clinical metadata used for grouping samples. Defaults to "Response".
                       #' @param stat_method String. The Method for comparing means.
                       #' @return A \code{ggplot} object.
                       plot_boxplot = function(dataset, feature_type, feature_name, group_col="Response",stat_method="wilcox.test") {

                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         feat_mat <- self$local_data[[dataset]]$Feature_Matrix[[feature_type]]
                         if (is.null(feat_mat)) stop("No matrix found for this feature type!")

                         if (!feature_name %in% feat_mat[,1]) stop(sprintf("Feature not found in matrix: %s", feature_name))
                         if (!group_col %in% colnames(clin_df)) stop(sprintf("Grouping column not found in clinical information: %s", group_col))

                         common_samples <- intersect(clin_df$Sample, colnames(feat_mat))
                         if (length(common_samples) == 0) stop("There is no matching sample ID between clinical information and feature matrix!")

                         plot_data <- data.frame(
                           Sample = common_samples,
                           Value = as.numeric(feat_mat[feat_mat[,1]==feature_name, common_samples]),
                           Group = as.character(clin_df[match(common_samples,clin_df$Sample), group_col])
                         )

                         plot_data <- plot_data[!is.na(plot_data$Group), ]

                         p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Value, fill = Group)) +
                           ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
                           ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
                           ggplot2::theme_classic() +
                           ggplot2::labs(title = sprintf("Feature: %s", feature_name),
                                         subtitle = sprintf("Dataset: %s", dataset),
                                         y = "Feature Value", x = group_col) +
                           ggplot2::theme(legend.position = "none")+
                           ggpubr::stat_compare_means(
                             method = stat_method,
                             label = "p.format",
                             label.x = 1.5
                           )

                         return(p)
                       },

                       #' @description Generate an ROC Curve to evaluate the predictive performance of a specific feature.
                       #' @param dataset String. The dataset identifier (e.g., "DTXR100001").
                       #' @param feature_type String. The type of feature matrix (e.g., "Gene Expression").
                       #' @param feature_name String. The specific feature name to evaluate (e.g., "CD274").
                       #' @param group_col String. The column name in the clinical metadata representing the true binary response. Defaults to "Response".
                       #' @param positive_class String. The label in the grouping column that represents the positive class (e.g., "R" for Responder). Defaults to "R".
                       #' @return A \code{ggplot} object.
                       plot_roc = function(dataset, feature_type, feature_name, group_col="Response", positive_class="R") {
                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         feat_mat <- self$local_data[[dataset]]$Feature_Matrix[[feature_type]]
                         if (is.null(feat_mat)) stop("No matrix found for this feature type!")

                         if (!feature_name %in% feat_mat[,1]) stop(sprintf("Feature not found in matrix: %s", feature_name))
                         if (!group_col %in% colnames(clin_df)) stop(sprintf("Grouping column not found in clinical information: %s", group_col))

                         common_samples <- intersect(clin_df$Sample, colnames(feat_mat))
                         if (length(common_samples) == 0) stop("There is no matching sample ID between clinical information and feature matrix!")

                         plot_data <- data.frame(
                           Sample = common_samples,
                           Value = as.numeric(feat_mat[feat_mat[,1]==feature_name, common_samples]),
                           Group = as.character(clin_df[match(common_samples,clin_df$Sample), group_col])
                         )
                         plot_data <- plot_data[!is.na(plot_data$Group), ]

                         if (!positive_class %in% plot_data$Group) stop("The specified 'positive_class' cannot be found!")
                         plot_data$Label <- ifelse(plot_data$Group == positive_class, 1, 0)

                         roc_obj <- pROC::roc(response = plot_data$Label, predictor = plot_data$Value, quiet = TRUE)
                         auc_val <- round(roc_obj$auc, 3)

                         p <- pROC::ggroc(roc_obj, color = "#d73027", size = 1) +
                           ggplot2::theme_minimal() +
                           ggplot2::geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey50") +
                           ggplot2::labs(title = sprintf("Feature: %s", feature_name),
                                         subtitle = sprintf("AUC = %s", auc_val),
                                         x = "Specificity", y = "Sensitivity")

                         return(p)
                       },

                       #' @description Generate a Pie Chart for categorical clinical information.
                       #' @param dataset String. The dataset identifier (e.g., "DTXR100001").
                       #' @param clinical_col String. The column name in the clinical metadata representing the categorical variable (e.g., "Response", "Sex").
                       #' @return A \code{ggplot} object.
                       plot_pie = function(dataset, clinical_col) {
                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         if (!clinical_col %in% colnames(clin_df)) stop(sprintf("Column not found in clinical information: %s", clinical_col))

                         clin_vector <- na.omit(clin_df[[clinical_col]])
                         if (length(clin_vector) == 0) stop("The specified clinical column only contains NA values!")

                         plot_data <- as.data.frame(table(clin_vector))
                         colnames(plot_data) <- c("Category", "Count")
                         plot_data$Percentage <- prop.table(plot_data$Count) * 100

                         plot_data$Label <- sprintf("%s\n(%.1f%%)", plot_data$Category, plot_data$Percentage)

                         p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = "", y = Count, fill = Category)) +
                           ggplot2::geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
                           ggplot2::coord_polar(theta = "y", start = 0) +
                           ggplot2::geom_text(ggplot2::aes(label = Label),
                                              position = ggplot2::position_stack(vjust = 0.5),
                                              size = 4, color = "white", fontface = "bold") +
                           ggplot2::theme_void() +
                           ggplot2::labs(title = sprintf("Distribution of %s", clinical_col),
                                         subtitle = sprintf("Dataset: %s (N = %d)", dataset, sum(plot_data$Count)),
                                         fill = clinical_col) +
                           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                                          plot.subtitle = ggplot2::element_text(hjust = 0.5))
                         return(p)
                       }
                     )
)
