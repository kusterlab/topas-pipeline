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

TARGET_CONFIG="config_growing.json"

while IFS= read -r config_path; do
  # Skip empty or comment lines
  [[ -z "$config_path" || "$config_path" =~ ^# ]] && continue

  echo "Processing config file: $config_path"

  if [[ ! -f "$config_path" ]]; then
    echo "Warning: Config file '$config_path' does not exist. Skipping."
    continue
  fi

  # Copy the config file to current directory with the fixed name
  cp "$config_path" "$TARGET_CONFIG"
  echo "Copied '$config_path' to '$TARGET_CONFIG'"

  # Run make (assuming make will pick up config_growing.json from current dir)
  echo "Running make command in background..."
  make CONFIG_FILE="$TARGET_CONFIG" MEMORY_LIMIT=200gb CPU_LIMIT=16 docker_all &
  PID=$!

  echo "Waiting for PID '$PID' to finish..."
  wait $PID

  echo "Finished run for: $config_path"
  echo "-----------------------------"
done < "$INPUT_FILE"
