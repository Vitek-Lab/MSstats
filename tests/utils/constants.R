#constants file 

#file paths of datasets on s3
datasets <- list(
    navarro_diaumpire_frags = "processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_FragSummary.RDS", 
    navarro_diaumpire_pepts = "processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_PeptideSummary.RDS", 
    navarro_diaumpire_prots = "processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_ProtSummary.RDS", 
    navarro_diaumpire_dia_annot = "processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_annotation.RDS", 
    navarro_openswath_input = "processed_data/processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_input.RDS", 
    navarro_openswath_annot = "processed_data/processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_annotation.RDS",
    datasets_rds = "datasets/datasets.RDS",
    metadata_rds = "datasets/metadata_new.RDS",
    processed_data = "processed_data"
)

results <- list(
    code_deploy_results = "code-deploy-results/"
)