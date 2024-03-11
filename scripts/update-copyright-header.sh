#!/bin/bash
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/update-copyright-header.sh
#
# Copyright (C) 2016-2024 Ian McCulloch <ian@qusim.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Research publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER

# This script automates the process of updating and managing file headers, including
# copyright notices, contributors' information, and license details.

# Lookup table for mapping contributor names to email addresses
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

# BuildBot's information.  We want to ignore commits from this user when scanning the logs
BUILDBOT_USER="MPToolkit BuildBot"
BUILDBOT_EMAIL="mptoolkit@qusim.net"

# Function to get the copyright strings for the header
function get_copyright_strings {
   local file=$1
   local original_header="$2"

   # Function to get a list of contributors for the copyright notice
   function get_contributors {
      local file=$1

      local contributors_git=$(git log --follow --invert-grep --author="$BUILDBOT_USER" --pretty=format:'%an' "$file" | sort -u) || exit 1
      local contributors_existing=$(echo "$original_header" | grep -oP "// Copyright \(C\) [0-9]+-[0-9]+ \K[^<]+(?= <)" | sort -u)

      echo -e "$contributors_git\n$contributors_existing" | grep -v '^$' | sort -u
   }

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

   # Function to get the dates of contributions for a specific contributor
   function get_contributor_dates {
      local file=$1
      local contributor="$2"

      # Get dates from existing copyright header
      local existing_dates=$(echo "$original_header" | grep -oP "// Copyright \(C\) [0-9]+(-[0-9]+)? $contributor" | grep -oP '[0-9]+(-[0-9]+)?' | sed 's/-/\n/g' | sort -u)

      # Get dates from git logs
      local git_dates=$(git log --follow --invert-grep --author="$BUILDBOT_USER" --pretty=format:'%an %ad' --date=format:%Y "$file" | awk -v contributor="$contributor" '$0 ~ contributor {print $NF}' | sort -u)

      # Combine and format dates from existing header and git logs
      local all_dates=$(echo -e "$existing_dates\n$git_dates" | grep -v '^$' | sort -u)
      local first_year=$(echo "$all_dates" | head -n 1)
      local last_year=$(echo "$all_dates" | tail -n 1)

      # Format the date range for the contributor
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
      copyright_strings+="Copyright (C) $(get_contributor_dates "$file" "$contributor") $contributor <$email>\n" || exit 1
   done <<< "$contributors"

   echo -ne "$copyright_strings"
}

# Function to retrieve the original header from a file
function get_original_header {

   # Check if the file exists
   if [ ! -f "$file" ]; then
      echo "Error: File not found: $file" >&2
      exit 1
   fi

   case "$file" in
   *.h | *.cpp | *.cc)
      local end_pattern="\/\/ ENDHEADER"
      ;;
   *.sh)
      local end_pattern="# ENDHEADER"
      ;;
   *.py)
      local end_pattern="# ENDHEADER"
      ;;
   *)
      echo "Unsupported file type: $file" >&2
      exit 1
      ;;
   esac

   # Use grep to check for the presence of ENDHEADER
   if grep -q "^$end_pattern" "$file"; then
      local original_header=$(awk "/$end_pattern/{print; exit} {print}" "$file")
      echo -e "$original_header"
   else
      echo ""
   fi
}

# Function to display usage information
function show_usage {
  cat <<EOF
Usage: $0 [--help|--show|--list|--diff|--commit|--backup|--exclude [path]|file [names...]]

Options:
  --help            : Display this help message.
  --show            : Display the updated headers.
  --list            : Don't update anything, just list the files to be examined.
  --add             : Add headers to files that don't currently have one.
  --diff            : Show the diff between the existing file content and the proposed changes.
  --apply           : Apply the changes to the headers.
  --commit          : Commit the changes to the git repository (implies --apply).
  --backup          : With --commit, save the old file as filename.bak.
  --exclude [path]  : Exclude this path from the search; can be used more than once.
  --file [names...] : Instead of searching, use this list of files.

EOF
}

# Function to generate the updated header template
function get_header {
   local file=$1
   local original_header=$2
   local top_directory=$(git rev-parse --show-toplevel)
   local relative_filename=$(realpath --relative-to="$top_directory" "$file")
   [[ -z "$relative_filename" ]] && exit 1

   # Check file extension to determine language
   case "$file" in
      *.h | *.cpp | *.cc)
      local comment="// "
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
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER"
      ;;
      *.sh)
      local comment="# "
      local header_template="#!/bin/bash
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# $relative_filename
#
# GENERATED_FILE_NOTICE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Research publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER"
      ;;
      *.py)
      local comment="# "
      local header_template="#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# $relative_filename
#
# GENERATED_FILE_NOTICE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Research publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER"
      ;;
      *)
      echo "Unsupported file type: $file" >&2
      exit 1
      ;;
   esac

   local copyright_strings=$(get_copyright_strings "$file" "$original_header" | sort) || exit 1
   [[ -z  "$copyright_strings" ]] && exit 1

   # Initialize the modified string
   local mcopyright_strings=""

   # Loop over each line and add the comment string
   while IFS= read -r line; do
      mcopyright_strings+="$comment$line\n"
   done <<< "$copyright_strings"

   copyright_strings=$(echo -e "$mcopyright_strings" | sed '$s/\\n$//')

   # Combine the updated copyright strings with the header template
   local updated_header=$(awk -v copyright="$copyright_strings" '/GENERATED_FILE_NOTICE/ { print copyright; next } 1' <<< "$header_template")


   #local updated_header=$(awk -v comment="$comment" -v copyright="$copyright_strings" '/GENERATED_FILE_NOTICE/ { print comment copyright; next } 1' <<< "$header_template")

   echo "$updated_header"
}

# Check if the current directory is a Git repository
if [ ! -d ".git" ]; then
    echo "Error: Not a Git repository. Exiting."
    show_usage
    exit 1
fi

# Parse command-line options
show=false
list=false
filelist=false
diff=false
apply=false
commit=false
backup=false
add=false
num_options=0
filenames=()
excludes=()

# Main loop to process files based on command-line options
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
    --list)
      ((num_options++))
      list=true
      ;;
    --backup)
      backup=true
      ;;
    --file)
      filelist=true
      ;;
    --diff)
      ((num_options++))
      diff=true
      ;;
    --apply)
      ((num_options++))
      apply=true
      ;;
    --commit)
      ((num_options++))
      commit=true
      if $apply; then
         ((num_options--))
      fi
      apply=true
      ;;
    --add)
      add=true
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
      if $filelist; then
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
   echo "Missing option, need one of --show|--list|--diff|--apply|--commit"
   show_usage
   exit 1
fi

if [ "$num_options" -ne 1 ]; then
  echo "Error: too many options specified."
  show_usage
  exit 1
fi

if $backup && ! $commit; then
   echo "Error: --backup option only makes sense with --apply."
   show_usage
   exit 1
fi

if $filelist; then
   find_command="printf '%s\n' \"\${filenames[@]}\""
else
   find_command="find . -type f \( -name '*.h' -o -name '*.cpp' -o -name '*.cc' -o -name '*.sh' -o -name '*.py' \) ${excludes[@]/#/-and -not -path }"
fi

# Check for uncommitted changes in the repository
if $apply && ([ -n "$(git diff --exit-code)" ] || [ -n "$(git diff --cached --exit-code)" ]); then
   read -p "There are uncommitted changes in the git repo. Do you want to continue? (y/N): " continue_confirm
   if [[ "$continue_confirm" != "y" ]]; then
      echo "Exiting without modifying headers."
      exit 1
   fi
fi

# track whether we actually changed a file using --commit
changedfiles=false

#Originally this was: eval "$find_command" | while IFS= read -r file; do
# but that runs in a subshell due to the | pipe, and we can't modify local variables inside a subshell
while IFS= read -r file; do

   # If --list is specified, just list the files
   if $list; then
      echo "$file"
      continue
   fi

   # scan the original header
   original_header=$(get_original_header "$file")

   # Generate the updated header template
   updated_header=$(get_header "$file" "$original_header") || exit 1

   # Show file content
   if $show; then
      echo -e "File: $file\n$updated_header\n"
      continue
   fi

   if [ "$updated_header" == "$original_header" ]; then
      echo "No changes required in $file"
      continue
   fi

   # Display the diff
   if $diff; then
      if [ -z "$original_header" ]; then
         echo "New header for file $file"
      fi
      diff -u -L "Current $file" -L "Updated $file" <(echo -e "$original_header") <(echo -e "$updated_header")
      echo
      continue
   fi

   if $apply; then

      if [ -z "$original_header" ] && ! $add; then
         echo "Ignoring $file as it doesn't have a header."
         continue
      fi

      changedfiles=true

      # Get the permissions of the original file
      original_mode=$(stat -c %a "$file")

      # Create a temporary file
      temp_file=$(mktemp) || exit 1

      # Add the updated header to the temp file
      echo -e "$updated_header" > "$temp_file"

      if [ -z "$original_header" ]; then
         if $add; then
            # if there is no header
            echo "Added new header to $file"

            # add a blank line, if there wasn't one already
            [[ -z $(head -n 1 "$file") ]] || echo >> "$temp_file"
            cat "$file" >> "$temp_file"
         else
            echo "oops, we shouldn't get here!" >&2
            exit 1
         fi
      else
         # there is an existing header
         echo "Updated header of $file"

         # Add the original file contents after // ENDHEADER or # ENDHEADER to the temp file
         awk '/^\/\/ ENDHEADER|^# ENDHEADER/{flag=1; next} flag' "$file" >> "$temp_file"
      fi

      # Backup the original file if --backup is specified
      if $backup; then
         mv "$file" "$file.bak"
      fi

      # Replace the old file with the updated one
      mv "$temp_file" "$file"
      # Update the new file's mode to match the original file
      chmod "$original_mode" "$file"
      if $backup; then
         echo "Old version saved as $file.bak"
      fi
   fi
done < <(eval "$find_command")

git_commit_command="git commit -a --author=\"$BUILDBOT_USER <$BUILDBOT_EMAIL>\" -m \"Auto-generate copyright headers\""

if $apply; then
   if $changedfiles || ([ -n "$(git diff --exit-code)" ] || [ -n "$(git diff --cached --exit-code)" ]); then
      if $commit; then
         # Commit the changes to the git repo
         echo -e "\nCommitting changes to the git repository..."
         eval "$git_commit_command"
         echo -e "\nSummary of the commit:"
         git log -1 --stat
      else
         echo -e "\nThere are changes to the headers. Commit these to the repository using:"
         echo -e "$git_commit_command"
         exit 1
      fi
   else
      echo -e "\nThere were no changes to the headers, nothing else to do."
      exit 1
   fi
fi
