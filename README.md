# Human eyeIntegration > v1.00 Shiny App
This is `pre-release` code and data. The datasets will undergo at least one more release.

# Web app (currently at 0.72)
[eyeIntegration.nei.nih.gov]()

# Install locally and run in three steps:
Warning! `get_eyeIntegration_datasets()` currently downloads ~9 GB of data and requires ~40 GB of space, uncompressed. 
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

# Where's the 40 GB of data going, by default?
```
system.file('app', package = 'eyeIntegrationApp')
```
