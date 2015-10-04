#!/bin/bash
# run this file after a compilation, it scans the .d files and produces a list of all of the source files that are used
cat *.d | sed 's/ \\$//' | tr ' ' '\n' | grep -v '^$' | grep -v ':' | sort | uniq

