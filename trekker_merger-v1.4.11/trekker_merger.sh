#!/bin/bash
set -e

#========================== Define Conda Paths =========================
PROFILE_CONDA_PATH=/opt/conda/envs/trekker
CONDA_SH_PATH=/opt/conda/etc/profile.d/conda.sh

#====================== Define Script Directory Paths ==================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#================== Function to validate samplesheet format ============
convert_to_linux_format() {
   local file=$1
   if grep -q $'\r$' "$file"; then
      echo "Samplesheet contains carriage return characters. Converting to Linux format..."
      sed -i 's/\r$//' "$file"
      echo "Success in removing carriage return characters. Samplesheet is updated to Linux format."
   else
      echo "Samplesheet is in Linux format."
   fi
   if [ -n "$(tail -c1 "$file")" ]; then
      echo >> "$file"
   fi
}

#======== Check if the required number of arguments is provided =========
if [ "$#" -ne 4 ]; then
   echo "Abort the trekker_merger. Error: Invalid input arguments"
   echo "Usage: bash $0 path_to_samplesheet.csv path_to_output_dir sample_ID profile"
   exit 1
fi

#================ Assign the input arguments to variables ================
SAMPLESHEET=$(realpath "$1")
OUT_DIR=$(realpath -m "$2")
SAMPLE_ID="$3"
PROFILE="$4"

#================ Check if samplesheet Exists & Extract Path =============
echo "Checking samplesheet..."
if [ ! -f "$SAMPLESHEET" ]; then
   echo "Error message:"
   echo "Abort the trekker_merger. Error: File '$SAMPLESHEET' does not exist. Please provide a valid path to samplesheet".
   exit 1
else
   SAMPLESHEET_DIR=$(dirname "$SAMPLESHEET")
   echo "'$SAMPLESHEET' exists"
   echo "Checking samplesheet format"
   convert_to_linux_format $SAMPLESHEET
fi

#================ Check if Output Dir Path Exists ================
echo "Checking Output directory path..."
BASE_OUT_PATH=$(dirname "$OUT_DIR")
if [ ! -d "$BASE_OUT_PATH" ]; then
   echo "Abort the trekker_merger. Error: The base output folder '$BASE_OUT_PATH' does not exist. Please provide a valid output directory path."
   exit 1
else
   echo "'$BASE_OUT_PATH exists'"
fi

#================ Check SAMPLE_ID ================
echo "Checking Output Sample ID prefix..."
if [[ "$SAMPLE_ID" =~ [^a-zA-Z0-9_] ]]; then
  SAMPLE_ID=$(echo "$SAMPLE_ID" | sed 's/[^a-zA-Z0-9_]/_/g')
  echo "Special characters found in $SAMPLE_ID. Using updated SAMPLE_ID: $SAMPLE_ID"
else
  echo "Using SAMPLE_ID: $SAMPLE_ID"
fi

#================ Check Profile and convert it to lowercase ================
PROFILE_UPDATED=$(echo "$PROFILE" | tr '[:upper:]' '[:lower:]')

#================ Define Directory Paths ================
INPUT_DIR=${SAMPLESHEET_DIR}/input/
LOG_DIR=${SAMPLESHEET_DIR}/log/${SAMPLE_ID}

#================ Check for existing Input Directory ================
if [ -d "$INPUT_DIR" ] && [ "$(ls -A "$INPUT_DIR")" ]; then
    echo "$INPUT_DIR exists. Please delete the folder or update its name and retrigger the trekker_merger."
    exit 1
fi

#================ Create Directories ================ 
mkdir -p "$INPUT_DIR"
mkdir -p "$OUT_DIR"
mkdir -p "${OUT_DIR}/intermediates/"
mkdir -p "$LOG_DIR"

#========================= Function to Handle Errors ========================
error_msg() {
   FAILED_COMMAND="$BASH_COMMAND"
   LOG_FILE="${FAILED_COMMAND##*/}"
   STEP_NAME="${LOG_FILE%.log}"
   echo "Abort the trekker_merger. Error: ${STEP_NAME} failed, Please refer to ${LOG_DIR}/${LOG_FILE} for details."
   exit 1
}

#======================= Extract samplesheet Columns =======================
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r id out_path; do
   seurat_file="$out_path/intermediates/${id}_seurat_spatial.rds"
   metrics_file="$out_path/${id}_summary_metrics.csv"
   if [ -f "$seurat_file" ] && [ -f "$metrics_file" ]; then   
      cp "$seurat_file" "$INPUT_DIR"
      cp "$metrics_file" "$INPUT_DIR"
      echo ${id}
      echo "seurat_file and metrics_file for ${id} are copied to a temporary input directory"
   else
      if [ ! -f "$seurat_file" ] || [ ! -f "$metrics_file" ]; then
         [ ! -f "$seurat_file" ] && echo "Seurat file not found for ${id}. Please check your samplesheet and provide a valid path to the output folder for sample ${id}."
         [ ! -f "$metrics_file" ] && echo "Metrics file not found for ${id}. Please check your samplesheet and provide a valid path to the output folder for sample ${id}"
         exit 1 
      fi
   fi	   
done

#======================= Fuction to check Tile_ID and Single_cell_assay =======================
check_tile_id_and_single_cell_consistency() {
   INPUT_DIR="$1"
   FIRST_FILE_CHECKED=false   
   for metrics_file in "$INPUT_DIR"/*_summary_metrics.csv; do
      if [ ! -f "$metrics_file" ]; then
         echo "Metrics file not found: $metrics_file"
         exit 1
      fi
      CURRENT_TILE_ID=$(grep 'Tile_ID' "$metrics_file" | cut -d',' -f2)
      CURRENT_SINGLE_CELL_ASSAY=$(grep 'Single_cell_assay' "$metrics_file" | cut -d',' -f2)
      if [ "$FIRST_FILE_CHECKED" = false ]; then
         TILE_ID="$CURRENT_TILE_ID"
	 SINGLE_CELL_ASSAY="$CURRENT_SINGLE_CELL_ASSAY"
	 FIRST_FILE_CHECKED=true
      else
         if [ "$TILE_ID" != "$CURRENT_TILE_ID" ] || [ "$SINGLE_CELL_ASSAY" != "$CURRENT_SINGLE_CELL_ASSAY" ]; then
            [ "$TILE_ID" != "$CURRENT_TILE_ID" ] && echo "Error: Trekker Tile ID must be the same for all samples. Please check your samplesheet and provide the 'output' folders for samples from the same Trekker tile."
            [ "$SINGLE_CELL_ASSAY" != "$CURRENT_SINGLE_CELL_ASSAY" ] && echo "Error: The Single_cell_assay must be the same for all samples. Please check your samplesheet and provide the 'output' folders for samples processed by the same single cell assay platform."
            exit 1
         fi
      fi 
   done
   if [ "$FIRST_FILE_CHECKED" = true ]; then
      echo "All samples in the samplesheet have the same Trekker Tile ID: "$TILE_ID" and are from the same Single_cell_assay platform: "$SINGLE_CELL_ASSAY""
   fi
}
check_tile_id_and_single_cell_consistency "$INPUT_DIR"

#======================= Running the trekker_merger pipeline  =======================
echo "Starting trekker_merger"

if [[ "$PROFILE_UPDATED" == "conda" ]]; then
   	source "$CONDA_SH_PATH"
   	conda activate "$PROFILE_CONDA_PATH"
	CONTAINER_CMD=""

elif [[ "$PROFILE_UPDATED" == "singularity" ]]; then
   	SINGULARITY_IMAGE=trekker-v1.4.11.sif
   	PROFILE_CONDA_PATH=/opt/conda/envs/trekker
	CONTAINER_CMD="singularity exec --bind ${SCRIPT_DIR},${INPUT_DIR},${OUT_DIR} ${SCRIPT_DIR}/${SINGULARITY_IMAGE}"

elif [[ "$PROFILE_UPDATED" == "docker" ]]; then
   	DOCKER_CONTAINER=tbusaspatial/trekker:v1.4.11
   	PROFILE_CONDA_PATH=/opt/conda/envs/trekker
	CONTAINER_CMD="docker run -v ${SCRIPT_DIR}:${SCRIPT_DIR} -v ${INPUT_DIR}:${INPUT_DIR} -v ${OUT_DIR}:${OUT_DIR} ${DOCKER_CONTAINER}"

else
   echo "Error message:"
   echo "Abort trekker_merger. Error: The 'profile' value argument is invalid. It must be one of the following: singularity, docker, or conda"
   exit 1
fi

echo "Running the trekker_merger using '$PROFILE'"

#=================== Set the trap to catch any ERR signal ===================
trap 'error_msg' ERR


echo "Start mergeSeurats"
${CONTAINER_CMD:+$CONTAINER_CMD }Rscript --vanilla ${SCRIPT_DIR}/mergeSeurats.R \
   ${INPUT_DIR} \
   spatial.rds \
   ${SAMPLE_ID} \
   ${OUT_DIR} \
   ${PROFILE_CONDA_PATH} >& ${LOG_DIR}/mergeSeurats.log
echo "mergeSeurats Done"

echo "Start mergeMetrics"
${CONTAINER_CMD:+$CONTAINER_CMD }python ${SCRIPT_DIR}/mergeMetrics.py \
   ${INPUT_DIR} \
   ${SAMPLE_ID} \
   ${OUT_DIR} >& ${LOG_DIR}/mergeMetrics.log
echo "mergeMetrics Done"

echo "Start genreport"
${CONTAINER_CMD:+$CONTAINER_CMD }Rscript --vanilla ${SCRIPT_DIR}/genreport.R \
   ${SAMPLE_ID} \
   ${TILE_ID} \
   ${SCRIPT_DIR}/htmlrender.Rmd \
   ${OUT_DIR} \
   ${SINGLE_CELL_ASSAY} \
   ${OUT_DIR}/${SAMPLE_ID}_summary_metrics_merged.csv \
   ${OUT_DIR}/${SAMPLE_ID}_variable_features_clusters_merged.csv \
   ${OUT_DIR}/${SAMPLE_ID}_variable_features_spatial_moransi_merged.txt \
   ${OUT_DIR}/intermediates/${SAMPLE_ID}_seurat_spatial_merged.rds \
   ${OUT_DIR}/${SAMPLE_ID}_ConfPositioned_seurat_spatial_merged.rds  >& ${LOG_DIR}/genreport.log
echo "genreport Done"
echo "Your trekker output has been merged. Please find the results in $OUT_DIR"

#======================= Remove intermediate input folder =======================
rm -r ${SAMPLESHEET_DIR}/input/

