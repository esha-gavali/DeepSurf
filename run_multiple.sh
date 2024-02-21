#!/bin/bash

# Define the directory containing your PDB files
PDB_DIR="receptor"

# Define the path to your model
MODEL_PATH="models"

# Define the base directory where you want to save output files
OUTPUT_BASE_DIR="output"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

start_time=$(date +%s)
# Loop through each PDB file in the directory
for PDB_FILE in "$PDB_DIR"/*.pdb; do
    # Extract the filename without the directory path and extension
    FILENAME=$(basename -- "$PDB_FILE" .pdb)
    
    # Define the output path for this specific PDB file
    OUTPUT_PATH="$OUTPUT_BASE_DIR/${FILENAME}_output"
    
    # Run the predict.py script with the current PDB file, model path, and output path
    python predict.py -p "$PDB_FILE" -mp "$MODEL_PATH" -o "$OUTPUT_PATH"
    
    # Optionally, echo a message indicating completion of this file
    echo "Completed: $PDB_FILE"
done

echo "All files processed."

end_time=$(date +%s)

# Calculate and print the total processing time
elapsed_time=$((end_time - start_time))
echo "Total processing time: $elapsed_time seconds."