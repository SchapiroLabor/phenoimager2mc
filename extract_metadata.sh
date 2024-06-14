#!/bin/bash

# Directory containing the TIFF files
input_dir="/Users/chiaraschiller/Documents/Schapiro/projects/SignalStar/data/Human_squamous_cell_carcinoma_stained_with_SignalStar_mIHC_technology/test"
# Output directory for the metadata files
output_dir="/Users/chiaraschiller/Documents/Schapiro/projects/SignalStar/data/Human_squamous_cell_carcinoma_stained_with_SignalStar_mIHC_technology/test/"

# Ensure the output directory exists
#mkdir -p "$output_dir"

# Loop through each TIFF file in the input directory
for tiff_file in "$input_dir"/*.tif; do
    # Get the base name of the TIFF file (without the directory path)
    base_name=$(basename "$tiff_file")
    # Create the output XML file name by replacing the .tif extension with .xml
    xml_file="$output_dir/${base_name%.tif}.xml"

    # Extract metadata and save it to the XML file
    identify -verbose "$tiff_file" > "$xml_file"

    # Print a message indicating that the metadata has been saved
    echo "Metadata for $tiff_file saved to $xml_file"
done