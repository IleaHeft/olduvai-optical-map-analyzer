1. Retrieve files from Pui Lab server
2. Upload files to Tesla
3. Make a folder with the sample name, move the zipped, tarball file into the folder, then unzip and extract the file
with ```gzip -d YOURFILENAME``` and ```tar -xf YOURFILENAME```
4. Process the structural variant outputs and run intersects with NBPF gene coordinates
  Inputs:
    * The folder that contains all of the data files (in folders by sample, unzipped, and extracted)
    * The coordinates (e.g. NBPF genes) that you want to run intersect on
    * The results folder
  Example usage -- remember to leave the "/" off the head of the data and results folders:
    ```bash scripts/analysis-of-duf-overlapping-sv.sh data/1000GenomesIrys/full-output-folders supporting-reference-data/duf1220-containing-gene-cords-combined-annotations.bed results/intersections-with-gene-coords/1000Genomes/2016-08-22```
5. Run an additional script that does some calculations on the results
  Inputs:
  * The file, generated from the script above, that gives the combined intersect results for all files
  * The results folder
  Example usage:
  ```bash scripts/perform-groupby-functions-on-combined-intersect-results.sh results/intersections-with-gene-coords/1000Genomes/2016-08-22/combined-intersect-results.bed results/intersections-with-gene-coords/1000Genomes/2016-08-22```
