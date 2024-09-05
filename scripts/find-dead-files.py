import os
import re


# List of directories to ignore (relative to the source directory)
IGNORE_DIRS = [
    "utils",
    "old",
    "attic",
    "misc",
    "junk",
    "linearalgebra/test",
    "quantumnumbers/test",
    "lattice/test",
    "models/old",
    "tensor/test",
    "pheap/test",
    "mpo/test",
    "common/test",
    "mp-algorithms/test",
    "linearalgebra/examples",
    "models/contrib/non-functional"
    # Add more directories as needed
]

# Function to parse the Makefile and extract the SRCDIR value
def parse_makefile_for_srcdir(makefile_path):
    with open(makefile_path, 'r') as makefile:
        for line in makefile:
            match = re.match(r'^SRCDIR\s*=\s*(.*)', line)
            if match:
                return match.group(1)
    return None

# Helper function to get the list of all .cpp files in the source directory, excluding ignored directories
def get_cpp_files(source_dir):
    cpp_files = []
    for root, _, files in os.walk(source_dir):
        # Skip any directories in the IGNORE_DIRS list
        relative_root = os.path.relpath(root, source_dir)
        if any(relative_root.startswith(ignore) for ignore in IGNORE_DIRS):
            continue
        for file in files:
            if file.endswith(".cpp"):
                cpp_files.append(os.path.join(root, file))
    return cpp_files

# Helper function to check if a corresponding .o file or executable exists in the build directory
def corresponding_object_or_executable_exists(cpp_file, source_dir, build_dir):
    # Get the base name of the file (without directory and extension)
    base_name = os.path.splitext(os.path.basename(cpp_file))[0]

    # Check for the object file
    object_file = base_name + ".o"
    object_exists = os.path.exists(os.path.join(build_dir, object_file))

    # Check for the executable (no extension)
    executable_file = base_name
    executable_exists = os.path.exists(os.path.join(build_dir, executable_file))

    # Return True if either the object file or executable exists
    return object_exists or executable_exists

# Main function to find dead cpp files
def find_dead_cpp_files(build_dir):
    makefile_path = os.path.join(build_dir, "Makefile")
    source_dir = parse_makefile_for_srcdir(makefile_path)

    if not source_dir:
        print("Source directory not found in Makefile.")
        return

    cpp_files = get_cpp_files(source_dir)
    dead_files = []

    for cpp_file in cpp_files:
        if not corresponding_object_or_executable_exists(cpp_file, source_dir, build_dir):
            dead_files.append(cpp_file)

    if dead_files:
        print("The following C++ files are not linked into any final programs:")
        for dead_file in dead_files:
            print(dead_file)
    else:
        print("No dead C++ files found.")

# Run the script
if __name__ == "__main__":
    build_dir = os.getcwd()  # Assumes the script is run from the build directory
    find_dead_cpp_files(build_dir)
