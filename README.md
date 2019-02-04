# Human eyeIntegration >= v1.00 Shiny App

# Web app
[https://eyeIntegration.nei.nih.gov](https://eyeIntegration.nei.nih.gov)

# Install locally and run in three steps:
Warning! `get_eyeIntegration_datasets()` currently downloads ~12 GB of data and requires ~45 GB of space, uncompressed. 
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
