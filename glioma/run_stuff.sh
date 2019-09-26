#!/bin/bash

#curated
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt curated_basis12.txt curated_results/ 2
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt curated_basis12.txt curated_results/ 5
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt curated_basis12.txt curated_results/ 10

#Immuno
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt immunoStates_basis.txt immunoStates_results/ 2
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt immunoStates_basis.txt immunoStates_results/ 5
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt immunoStates_basis.txt immunoStates_results/ 10

#LM22
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt LM22.txt lm22_results/ 2
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt LM22.txt lm22_results/ 5
Rscript glioma_ml_2.16.R genesymbols/Glioma_genesymbols.txt metadata/Glioma_metadata.txt LM22.txt lm22_results/ 10
