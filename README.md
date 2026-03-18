# DTXRquery

**DTXRquery** is an intuitive, R6-based client package designed to interact with the DisTxRESP database. It streamlines the entire bioinformatics workflow from data retrieval to visualization. 

With `DTXRquery`, you can easily:
* **Query & Filter**: Search the global metadata manifest using specific biological and clinical criteria (e.g., disease type, omics data, treatment regimens).
* **Batch Download**: Automatically download selected clinical data, feature matrices, and differential analysis results to your local machine.
* **In-Memory Analysis**: Load targeted datasets directly into the R environment as structured lists.
* **Built-in Visualization**: Generate high-quality, publication-ready plots (Volcano plots, Boxplots, and ROC curves) with simple, one-line commands.

## Installation

You can install the development version of DTXRquery from GitHub using the `devtools` or `remotes` package:

```R
# install.packages("devtools")
devtools::install_github("RuiyangZhai/DisTxRESP")
```

## Quick Start / Examples
Here is a typical workflow demonstrating how to initialize the client, filter data, load it into memory, and generate plots.

#### 1. Initialize the Client
Creating a new `DisTxRESP` object will automatically fetch the latest data manifest from the server.
```R
library(DTXRquery)

# Initialize the client
client <- DisTxRESP$new()

# View the full metadata
head(client$metadata)
```
#### 2. Filter and Load Data
You can chain methods together to filter the datasets you need and load them directly into memory for analysis.
```R
# Filter specific data sets and feature types
client$filter_metadata(
  dtxr = c("DTXR100001","DTXR100002"),
  feature = c("Clinical Data","Pathway Activity","Gene Expression")
)

# Check the filtered subset table
head(client$sub_table)

# Batch download files
client$download_files(output_dir = "my_data")

# Load into memory
client$load_to_memory(local_dir = "my_data")
```
_(Note: If you do not want to save the file locally, you can also directly use `client$load_to_memory()`.)_
#### 3. Data Visualization
Once the data is loaded into memory, you can use the built-in visualization methods.

__Volcano Plot__

Visualize differential expression results quickly:
```R
# Plot a volcano plot for a specific dataset
volcano_plt <- client$plot_volcano(
  dataset = "DTXR100002", 
  feature_type = "Gene Expression",
  logFC_col = "LogOR",
  pval_col = "P.value",
  fc_cutoff = 1.0, 
  p_cutoff = 0.05
)
print(volcano_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/CIsubtype.png?raw=true" width="600">
</div>

__Boxplot__

Compare the expression or values of a specific feature across clinical groups:
```R
# Plot the expression of CD274 (PD-L1) across treatment response groups
box_plt <- client$plot_boxplot(
  dataset = "DTXR100002", 
  feature_type = "Gene Expression", 
  feature_name = "CD274", 
  group_col = "Response"
)
print(box_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/CIsubtype.png?raw=true" width="600">
</div>

__ROC Curve__

Evaluate the predictive performance of a biomarker:
```R
# Generate an ROC curve to assess CD274 as a predictor for Responder ("R")
roc_plt <- client$plot_roc(
  dataset = "DTXR100001", 
  feature_type = "Gene Expression", 
  feature_name = "CD274", 
  group_col = "Response", 
  positive_class = "R"
)
print(roc_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/CIsubtype.png?raw=true" width="600">
</div>

## Contact
Any technical question please contact Ruiyang Zhai (zhairuiyang@foxmail.com) or Te Ma (mate.compbio@foxmail.com).
