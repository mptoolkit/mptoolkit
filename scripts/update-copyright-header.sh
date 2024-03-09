#!/bin/bash

# Lookup table for names to email addresses
declare -A email_lookup=(
  ["Ian McCulloch"]="ian@qusim.net"
  ["Jesse Osborne"]="j.osborne@uqconnect.edu.au"
  ["Seyed Saadatmand"]="s.saadatmand@uq.edu.au"
  ["Tomohiro"]="tomohiro.hashizume@uq.net.au"
  ["Henry Nourse"]="henry.nourse@uqconnect.edu.au"
  ["Stefan Depenbrock"]="Stefan.Depenbrock@physik.uni-muenchen.de"
  ["Fei Zhan"]="enfeizhan@gmail.com"
  # Add more entries as needed
)

# Function to get the email address for a contributor
function get_email {
  local contributor="$1"

  if [ "${email_lookup[$contributor]+isset}" ]; then
    echo "${email_lookup[$contributor]}"
  else
    echo "Error: Email address for $contributor not found in email lookup table." >&2
    echo "Names appearing in the git log are: " >&2
    git log --all --format="%aN <%aE>" | grep "$contributor" | sort -u >&2
    exit 1
  fi
}

# Function to get the copyright strings for the header
function get_copyright_strings {
   local file=$1

   # Function to get a list of contributors for the copyright notice
   function get_contributors {
      local file=$1

      local contributors_git=$(git log --follow --pretty=format:'%an' "$file" | sort -u) || exit 1
      local contributors_existing=$(grep -oP "// Copyright \(C\) [0-9]+-[0-9]+ \K[^<]+|// Copyright \(C\) [0-9]+-[0-9]+ \K[^<]+(?= <)" "$file" | sed 's/ $//' | sort -u) || exit 1

      echo -e "$contributors_git\n$contributors_existing" | grep -v '^$' | sort -u
   }

   # Function to get the dates of contributions for a specific contributor
   function get_contributor_dates {
      local file=$1
      local contributor="$2"

      # Get dates from existing copyright header
      local existing_dates=$(grep -oP "// Copyright \(C\) [0-9]+(-[0-9]+)? $contributor" "$file" | grep -oP '[0-9]+(-[0-9]+)?' | sed 's/-/\n/g' | sort -u)

      # Get dates from git logs
      local git_dates=$(git log --follow --pretty=format:'%an %ad' --date=format:%Y "$file" | grep "$contributor" | awk '{print $NF}' | cut -d'-' -f1 | sort -u)

      # Combine and format dates
      local all_dates=$(echo -e "$existing_dates\n$git_dates" | grep -v '^$' | sort -u)
      local first_year=$(echo "$all_dates" | head -n 1)
      local last_year=$(echo "$all_dates" | tail -n 1)

      if [ "$first_year" = "$last_year" ]; then
         echo "$first_year"
      else
         echo "$first_year-$last_year"
      fi
   }

   local contributors=$(get_contributors "$file") || exit 1
   local copyright_strings=""

   while IFS= read -r contributor; do
      local email=$(get_email "$contributor") || exit 1
      if [ -z "$email" ]; then
         echo "Error occurred while processing file $file" >&2
         exit 1
      fi
      copyright_strings+="// Copyright (C) $(get_contributor_dates "$file" "$contributor") $contributor <$email>\n" || exit 1
   done <<< "$contributors"

   echo -e "$copyright_strings"
}

# Function to get the original header from the file
function get_original_header {
  local file="$1"

  # Check if the file exists
  if [ ! -f "$file" ]; then
    echo "Error: File not found: $file" >&2
    exit 1
  fi

  # Get the header content up to // ENDHEADER (inclusive)
  original_header=$(awk '/\/\/ ENDHEADER/{print; exit} 1' "$file")

  # If there is no header, exit an empty string
  if [ -z "$original_header" ]; then
    echo ""
  else
    echo "$original_header"
  fi
}

# Function to display usage information
function show_usage {
  cat <<EOF
Usage: $0 [--help|--show|--files|--exclude [path]|--diff|--commit|file [names...]]

Options:
  --help            : Display this help message.
  --show            : Display the updated header for each file (default action with --files).
  --files           : List the files to be examined.
  --exclude [path]  : Exclude this path from the search; can be used more than once.
  --diff            : Show the diff between the existing file content and the proposed changes.
  --commit          : Apply the changes to the header.
  --backup          : With --commit, save the old file as filename.bak.
  --file [names...] : Show the updated header for this list of files. Can be used with --diff.

EOF
}

# Function to get the first and last years of contribution for each contributor
function get_header {
  local file=$1
  local top_directory=$(git rev-parse --show-toplevel)
  local relative_filename=$(realpath --relative-to="$top_directory" "$file")
  [[ -z  "$relative_filename" ]] && exit 1
  local header_template="// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// $relative_filename
//
// GENERATED_FILE_NOTICE
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER"

   local copyright_strings=$(get_copyright_strings "$file") || exit 1
   [[ -z  "$copyright_strings" ]] && exit 1
   local updated_header=$(awk -v copyright="$copyright_strings" '/GENERATED_FILE_NOTICE/ { print copyright; next } 1' <<< "$header_template")

   echo "$updated_header"
}

# Parse command-line options
show=false
files=false
showfile=false
diff=false
commit=false
backup=false
num_options=0
filenames=()
excludes=()

while [ "$#" -gt 0 ]; do
  case "$1" in
    --help)
      show_usage
      exit 0
      ;;
    --show)
      ((num_options++))
      show=true
      ;;
    --files)
      ((num_options++))
      files=true
      ;;
    --backup)
      backup=true
      ;;
    --file)
      #((num_options++))
      showfile=true
      ;;
    --diff)
      ((num_options++))
      diff=true
      ;;
    --commit)
      ((num_options++))
      commit=true
      ;;
    --exclude)
      if [ "$#" -ge 2 ]; then
         excludes+=("'./$2*'")
         shift
      else
         echo "Error: --exclude option requires a directory or file argument." >&2
         show_usage
         exit 1
      fi
      ;;
    *)
      if $showfile; then
        filenames+=("$1")
      else
        show_usage
        exit 1
      fi
      ;;
  esac
  shift
done

if [ "$num_options" -eq 0 ]; then
   # If we've specified some files, then default to --show
   if $showfile; then
      show=true
      num_options=1
   else
      show_usage
      exit 1
   fi
fi

if [ "$num_options" -ne 1 ]; then
  echo "Error: too many options specified."
  show_usage
  exit 1
fi

if $backup && ! $commit; then
   echo "Error: --backup option only makes sense with --commit."
   show_usage
   exit 1
fi

if $showfile; then
   find_command="printf '%s\n' \"\${filenames[@]}\""
else
   find_command="find . -type f \( -name '*.h' -o -name '*.cpp' -o -name '*.cc' \) ${excludes[@]/#/-and -not -path }"
fi

eval "$find_command" | while IFS= read -r file; do

  # If --files is specified, just list the files
  if $files; then
    echo "$file"
    continue
  fi

  # Check if the file contains the marker "// ENDHEADER"
  if grep -q '// ENDHEADER' "$file"; then
    # Generate the updated header template
    updated_header=$(get_header "$file") || exit 1

    # Show file content
    if $show; then
      echo -e "File: $file\n$updated_header\n"
    fi

    # Display the diff
    if $diff; then
      original_header=$(get_original_header "$file") || exit 1
      echo -e "File: $file\nDifferences in the updated header:"
      diff -u -L "Original $file" -L "Updated $file" <(echo -e "$original_header") <(echo -e "$updated_header")
      echo
    fi

    # Apply changes only if --commit is specified
    if $commit; then
      # Create a temporary file
      temp_file=$(mktemp) || exit 1


      original_header=$(get_original_header "$file")
      if [ "$updated_header" == "$original_header" ]; then
         echo "No changes required in $file"
         continue
      fi

      # Add the updated header to the temp file
      echo -e "$updated_header" > "$temp_file"

      # Add the original file contents after // ENDHEADER to the temp file
      sed -n '/\/\/ ENDHEADER/,$p' "$file" | tail -n +2 >> "$temp_file"

      # Backup the original file if --backup is specified
      if $backup; then
         mv "$file" "$file.bak"
      fi

      # Move the temp file to the actual filename
      mv "$temp_file" "$file"
      echo "Updated header in $file"
      if $backup; then
         echo "Old version saved as $file.bak"
      fi
    fi
  else
    # Add new header only if --commit is specified
    if $commit; then
      # Create a temporary file to write the updated content
      temp_file=$(mktemp) || exit 1

      # Add new header to the temporary file
      new_header=$(get_header "$file") || exit 1
      echo -e "$new_header\n$(cat "$file")" > "$temp_file"

      # Backup the original file if --backup is specified
      if $backup; then
        mv "$file" "$file.bak"
      fi

      # Replace the old file with the updated one
      mv "$temp_file" "$file"
      echo "Added new header to $file"
      if $backup; then
         echo "Old version saved as $file.bak"
      fi
    fi
  fi
done
