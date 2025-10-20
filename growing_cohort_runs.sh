#!/bin/bash

# growing_cohort_runs.sh
# Usage: ./growing_cohort_runs.sh config_input_file.txt

INPUT_FILE="$1"

if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 config_input_file.txt"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: File '$INPUT_FILE' not found."
  exit 2
fi

DELIM=$'\t'  # assume tab-delimited

# Find index of the 'config_file' column
#HEADER=$(head -n 1 "$INPUT_FILE")
#IFS="$DELIM" read -ra COLS <<< "$HEADER"

#CONFIG_INDEX=-1
#for i in "${!COLS[@]}"; do
#  if [[ "${COLS[$i]}" == "config_file" ]]; then
#    CONFIG_INDEX=$i
#    break
#  fi
#done

#if [[ $CONFIG_INDEX -eq -1 ]]; then
#  echo "Error: Column 'config_file' not found in header."
#  exit 3
#fi



# Process each line starting from line 2
tail -n +2 "$INPUT_FILE" | while IFS="$DELIM" read -ra LINE; do
  CONFIG_PATH="${CONFIG_PATH%/}"

  CONFIG_NAME="$(basename "$CONFIG_PATH")"
  #echo CONFIG_NAME
  DEST_PATH="$SCRIPT_DIR/$CONFIG_NAME"
  echo "$CONFIG_NAME"
  echo "$CONFIG_PATH"
  #echo "Copying $CONFIG_PATH to $DEST_PATH"
  #cp "$CONFIG_PATH" "$DEST_PATH"
  echo "$DEST_PATH"
  

  echo "Running: CONFIG_FILE=$DEST_PATH MEMORY_LIMIT=200gb CPU_LIMIT=16 make docker_all"
  #CONFIG_FILE="$DEST_PATH" MEMORY_LIMIT=200gb CPU_LIMIT=16 make docker_all
done
