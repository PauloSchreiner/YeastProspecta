import os
import re
import argparse

# To test, run this from the snakemake root dir:
# python3 workflow/scripts/utils/rename_samples.py --dir data/raw_reads

# To execute:
# python3 workflow/scripts/utils/rename_samples.py --dir data/raw_reads --execute


# ==========================================
# REGEX EXPLANATION
# ==========================================
# Pattern breakdown for: 'SF56-2 NL4_F06.ab1'
# ^([A-Za-z]+) : Group 1 (Prefix) -> SF
# (\d+)        : Group 2 (Number) -> 56
# (-\d+)?      : Group 3 (Suffix) -> -2 (Optional)
# \s+          : Matches the space(s) between the sample name and the well code
# (NL1|NL4)    : Group 4 (Direction) -> NL1 or NL4
# _[A-Za-z]\d+ : Matches the well code (e.g., _F06) to discard it
# \.(ab1|seq)$ : Group 5 (Extension) -> .ab1 or .seq

REGEX_PATTERN = r"^([A-Za-z]+)(\d+)(-\d+)?\s+(NL1|NL4)_[A-Za-z]\d+\.(ab1|seq)$"


def rename_files(directory: str, dry_run: bool = True):
    """
    Iterates through the target directory, matches sequencing files against 
    the expected pattern, standardizes their names (z-fill and direction), 
    and renames them.
    """
    
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' not found.")
        return

    processed_count = 0
    skipped_count = 0

    print(f"--- Starting File Renaming Process ---")
    print(f"Directory: {directory}")
    print(f"Mode: {'DRY-RUN (No files will be changed)' if dry_run else 'ACTIVE (Files WILL be renamed)'}\n")

    for filename in os.listdir(directory):
        if not filename.endswith(('.ab1', '.seq')):
            continue

        match = re.search(REGEX_PATTERN, filename)

        if match:
            # Extract groups
            prefix = match.group(1)            
            number = match.group(2)            
            suffix = match.group(3) or ""      
            direction_raw = match.group(4)     
            extension = match.group(5)         

            # --- TRANSFORMATIONS ---
            number_padded = number.zfill(3)

            if direction_raw == "NL1":
                direction_clean = "F"
            elif direction_raw == "NL4":
                direction_clean = "R"
            else:
                direction_clean = direction_raw 

            new_filename = f"{prefix}{number_padded}{suffix}_{direction_clean}.{extension}"

            # --- EXECUTION ---
            old_filepath = os.path.join(directory, filename)
            new_filepath = os.path.join(directory, new_filename)

            if dry_run:
                print(f"[DRY-RUN] '{filename}'  --->  '{new_filename}'")
            else:
                os.rename(old_filepath, new_filepath)
                print(f"[RENAMED] '{filename}'  --->  '{new_filename}'")
            
            processed_count += 1
            
        else:
            print(f"[SKIPPED] '{filename}' does not match the expected pattern.")
            skipped_count += 1

    print(f"\n--- Summary ---")
    print(f"Files correctly formatted: {processed_count}")
    print(f"Files skipped (unmatched): {skipped_count}")
    
    if dry_run:
        print("\n[NOTE] This was a DRY-RUN. To apply changes, run the script with the --execute flag.")


if __name__ == "__main__":
    # Setup argparse for terminal execution
    parser = argparse.ArgumentParser(
        description="Standardizes sequencing raw files (.ab1, .seq) by zero-filling numbers and converting NL1/NL4 to F/R."
    )
    
    # Argument 1: The directory (Required)
    parser.add_argument(
        "-d", "--dir", 
        type=str, 
        required=True, 
        help="Path to the directory containing the raw sequencing files."
    )
    
    # Argument 2: The execute flag (Optional)
    parser.add_argument(
        "--execute", 
        action="store_true", 
        help="Actually rename the files on disk. If omitted, performs a safe DRY-RUN."
    )
    
    # Parse the arguments typed in the terminal
    args = parser.parse_args()
    
    # If --execute is present in the terminal, args.execute is True (so dry_run becomes False)
    is_dry_run = not args.execute
    
    # Call the main function
    rename_files(directory=args.dir, dry_run=is_dry_run)