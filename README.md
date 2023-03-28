# Human eyeIntegration >= v2.0 Shiny App

# Updated to serve the 2022 data update. New samples. Metadata cleaning. Recount3-based quantification

# Web app
[https://eyeIntegration.nei.nih.gov](https://eyeIntegration.nei.nih.gov)

# scEiaD based sqlite file:
https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/scEiaD.2023_03_02.sqlite.gz

# Files for use in PCA Visualization can be found here:

https://hpc.nih.gov/~parikhpp/eyeIntegration_app_v2/eye_pca_data.RData.gz
https://hpc.nih.gov/~parikhpp/eyeIntegration_app_v2/eye_pca_data_with_GTEx.RData.gz
https://hpc.nih.gov/~parikhpp/eyeIntegration_app_v2/eye_percentVar_data.RData.gz

Updated PCA Files (03-28-23)
https://hpc.nih.gov/~parikhpp/eyeIntegration_app_v2/EiaD_pca_analysis.Rdata.gz
https://hpc.nih.gov/~parikhpp/eyeIntegration_app_v2/eyeIntegration_2023_pca.RData.gz

# All below is defunct right now

# Install locally and run in three steps:
Warning! `get_eyeIntegration_datasets()` currently downloads ~20 GB of data and requires ~65 GB of space, uncompressed. 
```
devtools::install_github('davemcg/eyeIntegration_app')
eyeIntegrationApp::get_eyeIntegration_datasets()
eyeIntegrationApp::run_eyeIntegration()
```

# Install locally (custom location!) and run in three steps:
```
devtools::install_github('davemcg/eyeIntegration_app')
eyeIntegrationApp::get_eyeIntegration_datasets(destdir='/your/path/to/app')
eyeIntegrationApp::run_eyeIntegration(app_path='/your/path/to/app')
```

# Already have the most recent SQLite databases and just need new web app resources?
```
devtools::install_github('davemcg/eyeIntegration_app')
eyeIntegrationApp::get_eyeIntegration_datasets(file = "eyeIntegration_v101_01_lightweight.tar.gz",
                                                       destdir='/your/path/to/app')
eyeIntegrationApp::run_eyeIntegration(app_path='/your/path/to/app')
```

# Where's the 45 GB of data going, by default?
```
system.file('app', package = 'eyeIntegrationApp')
```
