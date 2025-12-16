#!/bin/bash

#SBATCH --job-name tarball
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=56GB
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/tarball/SLURM-%j.out
#SBATCH --error=slurm-logs/tarball/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Get arguments
DIR_TO_COMPRESS=$1
NEW_NAME=$2

# Check arguments provided
if [ -z "$DIR_TO_COMPRESS" ] || [ -z "$NEW_NAME" ]; then
    echo "Error: Missing arguments"
    echo "Usage: sbatch tarball.sh <directory> <output_path_name.tar.gz>"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$NEW_NAME")
mkdir -p "$OUTPUT_DIR"

echo "Compressing: $DIR_TO_COMPRESS"
echo "Output: $NEW_NAME"
echo "Using 8 threads with pigz"

# Compress with pigz using 8 threads
tar --use-compress-program='pigz -p 16 --best' -cf "$NEW_NAME" "$DIR_TO_COMPRESS"

echo "Compression complete!"
echo "Output file: $NEW_NAME"
ls -lh "$DIR_TO_COMPRESS"
ls -lh "$NEW_NAME"