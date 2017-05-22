# irys-duf1220

irys-duf1220 is a set of tools for utilizing data produced by the Irys system to analyze structural variation in DUF1220 containing genes

# Getting Started
You will want to run 2 types of analysis.  

You will want to analyze the results of **molecule to reference** alignments, with multimatch on.  

You will also want to analyze the **molecule to contig - contig to reference** alignments, in order to recover information for certain genes (e.g. NBPF19) for which very little information is recovered with the molecule to reference approach.  There are two steps, the first is to analyze the contig to reference data.  The second is to analyze the molecule to contig data.  The results of both of these steps are files that are equivalent to the "con2-con3-distance" files produced by the "standard", molecule to reference analysis.

# Analyze the data
All of the necessary code has been packaged into one script, _irys-analyzer.sh_ with an accompanying _config.sh_ script to set various parameters and file paths.  The script will automatically analyze both HLS region and CON1 region variation. Depending on the **alignment_type** you set in the _config.sh_ script, the _irys-analyzer.sh_ script will run the appropriate analysis.  


For a "complete" analysis of a set of samples, you will want to run the _irys-analyzer.sh_ three separate times, each time specifiying a different alignment type in the _config.sh_ script.  

When you swith between running the MolRef and the MolContig or ContigRef data, make sure you also update the data folder (the MolRef data is likely in one folder and the MolContig and ContigRef data are in a separate folder).    

1) Ensure parameters in the config.sh file are set as desired
- alignment type (one of: MolRef, MolContig, or ContigRef)
- sample directory
- file name setup
- number of samples
- conf_spread
  - the **conf_spread** parameter defines the minimum difference between the maximum confidence score and other confidence scores for the molecule for the alignment to be retained
- link_dist
  - the **link_dist** parameter controls the size of the "bin" in the peak-caller.sh script 
- min_mols_in_cluster
  - the **min_mols_in_cluster** parameter defines how many molecules you want in each cluster in order to call the cluster as an SV 

2) Run the script:
```
bash irys-analyzer.sh
```
**Notes:**
The script automatically generate two sets of results, one in which molecules with a nick between the aligned con2 and con3 nicks are filtered out (the "adj-mols-only" folder within the main results folder), and a set of results in which that filter is not applied (the main results folder)  

File structure:
- For MolRef data, the script is set up to operate on a single folder that contains all of the files (xmap, r.cmap, and q.cmap) for every sample





