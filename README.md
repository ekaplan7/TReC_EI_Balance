# TReC_EI_Balance

# Overview
Comprehensive Behavioral Intervention for Tics (CBIT) is the gold standard treatment for Tourette Syndrome (TS), however roughly half of subjects don't show treatment response. Recent work suggests that the aperiodic exponent of the EEG power spectrum may serve as a non-invasive proxy for cortical excitation/inhibtion (E/I) balance.  Flatter slopes hae been associated with relatively greater excitation, whereas steeper slopes are thought to reflect relatively greater inhibition. Individuals with TS show heightened aperiodic activity during sensorimotor processing. This project investigates whether CBIT alters E/I balance in the brains of adolescents with TS. 

# Conceptual Framework
This project is grounded in three core ideas: 
- TS has been linked to altered inhibitory control and increased neural noise within cortico-striatal-thalamo-cortical circuits
- The aperiodic exponent has been proposed as a non-invasive marker related to E/I balance
- CBIT trains behavioral inhibtion, and therefore may influence underlying cortical excitability

The primary objective is to examing whether aperiodic parameters shift following treatment and whether these shifts relate to tic severity.

# Repository Layout
  TRec_EI_Preproc.m -> original preprocessing script. Takes in both tic onset and offset codes. 
        # File is organized such that you can use command + F to search for 'UPDATE' which will show all places in the script that must be updated for subj-specific data . 

  TReC_EI_Preproc_Final.m -> Current working preprocessing code. Takes in ONLY onset markers. Offset markers are estimated based on amplitude thresholding. 
        # File is similarly organized, you can search for 'UPDATE' to see blocks where subj-specific data is required 

  ts-ei-analysis.ipynb -> SpecParam workflow. Takes in folder of preprocessed subjects PRE and POST, runs specparam and plots results: individual and grand average PSDs and fitted model, trajectories, etc 

# Data 
Data for this project are collected as part of the broader TReC study in collaboration between UCSD and the Univeristy of Minnesota (UMN). This multi-site design enables larger and more diverse samples. 

Participants completed assessments both before and after an 8-week CBIT intervention. 

EEG data are organized by participant and session, with associated metadata including demographic ariables and clinical severity scores. fMRI and behavioral data are maintained in the larger TReC study and may be integrated for multimodal analyses in the future. 

# Data Access & Privacy
The data used in this project are derived from a clinical research study involving pediatric participants and therefore contain protected health information (PHI). For this reason, raw and processed data are not publicly available through this repository. All data collection and handling procedures are conducted in accordance with institutional review board (IRB) approvals at both collaborating sites and comply with applicable data protection and privacy regulations. Researchers interested in accessing the dataset for collaboration or secondary analyses may contact the study team to inquire about data availability. Any data sharing would require appropriate institutional approvals, data use agreements, and compliance with ethical and regulatory standards.
