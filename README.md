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
devtools::install_github("RuiyangZhai/DTXRquery")
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
You can chain methods together to filter the required data sets, download them and load them into memory for analysis.
```R
# Filter specific data sets and feature types
client$filter_metadata(
  dtxr = c("DTXR600031","DTXR500068"),
  feature = c("Clinical Data","Microbial Abundance")
)

# Check the filtered subset table
head(client$sub_table)

# Batch download files
client$download_files(output_dir = "my_data")

# Load into memory
client$load_to_memory(local_dir = "my_data")
```

#### 3. Data Visualization
Once the data is loaded into memory, you can use the built-in visualization methods.

__Pie Plot__
Pie chart for visualizing clinical baseline characteristics：
```R
# View the proportion of different response groups in the dataset
pie_plt <- client$plot_pie(dataset = "DTXR600031",clinical_col = "Response")

# View gender ratio of patients
# pie_plt <- client$plot_pie(dataset = "DTXR600031",clinical_col = "Sex")
print(pie_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/pie_plt.png?raw=true" width="600">
</div>

__Volcano Plot__

Visualize differential expression results quickly:
```R
# Plot a volcano plot for a specific dataset
volcano_plt <- client$plot_volcano(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance",
  logX_col = "LogOR",
  pval_col = "P.value",
  fc_cutoff = 1.0, 
  p_cutoff = 0.05
)
print(volcano_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/volcano_plt.png?raw=true" width="600">
</div>

__Boxplot__

Compare the expression or values of a specific feature across clinical groups:
```R
# Plot the abundance of Akkermansia across treatment response groups
box_plt <- client$plot_boxplot(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance", 
  feature_name = "g__Akkermansia", 
  group_col = "Response"
)
print(box_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/box_plt.png?raw=true" width="600">
</div>

__ROC Curve__

Evaluate the predictive performance of a biomarker:
```R
# Generate an ROC curve to assess Akkermansia as a predictor for Responder ("R")
roc_plt <- client$plot_roc(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance", 
  feature_name = "g__Akkermansia", 
  group_col = "Response", 
  positive_class = "R"
)
print(roc_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/roc_plt.png?raw=true" width="600">
</div>

## Contact
Any technical question please contact Ruiyang Zhai (zhairuiyang@foxmail.com) or Te Ma (mate.compbio@foxmail.com).
