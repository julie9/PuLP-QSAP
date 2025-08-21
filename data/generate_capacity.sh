#!/bin/bash

# Usage:
#    ./generate_capacity.sh 4 5 4 1
# This will create a file named processor_capacity_4x5_4x1.txt with the contents:
# 5
# 5
# 5
# 5
# 1
# 1
# 1
# 1

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ] || [ $(($# % 2)) -ne 0 ]; then
  echo "Usage: $0 <repetitions1> <number1> [<repetitions2> <number2> ...]"
  exit 1
fi

# Create a single temporary file
temp_file=$(mktemp)

# Initialize a summary string for the filename
summary=""

# Iterate over the arguments in pairs
while [ "$#" -gt 0 ]; do
  repetitions=$1
  number=$2
  shift 2

  # Append the number to the temp file the specified number of times
  for ((i = 0; i < repetitions; i++)); do
    echo "$number" >> "$temp_file"
  done

  # Append to the summary string
  summary+="${repetitions}x${number}_"
done

# Remove the trailing underscore from the summary
summary=${summary%_}

# Rename the temp file to include the summary in the filename
output_file="processor_capacity_${summary}.txt"
mv "$temp_file" "$output_file"

echo "Numbers written to $output_file"