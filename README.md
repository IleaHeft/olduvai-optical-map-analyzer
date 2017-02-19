# irys-duf1220

irys-duf1220 is a set of tools for utilizing data produced by the Irys system to analyze structural variation in DUF1220 containing genes

# Getting Started

All of the individual scripts have been packaged into a single script, _irys-analyzer.sh_ with an accompanying _config.sh_ script to set various parameters and file paths.  The script will automatically analyze both HLS region and CON1 region variation.

1) Ensure parameters in the config.sh file are set as desired  
>Note that the **link_dist** parameter controls the size of the "bin" in the peak-caller.sh script, and changing this parameter could change your results.  

2) Run the script:
```
bash irys-analyzer.sh
```

# File name formatting  
If the file names use underscore separators, convert them to dashes with _renamefiles.sh_.  First move the script into the folder with the files that need the name change,then cd into that directory and run as shown below.  Then move the script back to the scripts directory. Note: if you later intend to open the xmap with the changed file name, you'll need to revert the file names of the r.cmap and q.cmap files to match their original formatting - if you do not, then IrysView will throw an error.
  
 ```
 bash renamefiles.sh
 ```

