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

# Where's the 40 GB of data going?
```
system.file('app', package = 'eyeIntegrationApp')
```
