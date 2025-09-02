# -*- coding: utf-8 -*-
import zipfile
import os
from io import StringIO
import argparse
import copy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import shutil


def determine_type(value):
    # Try to cast to an int
    try:
        int(value)
        return "int"
    except ValueError:
        # Not an int
        pass
    return "string"

def is_fasta_file(input_arg):
    """Basic check to see if a file is a FASTA file based on its first line."""
    if input_arg.endswith('.fasta') or input_arg.endswith('.fa') or input_arg.endswith('.fna'):
        return True
    else:
        return False

def is_fastq_file(input_arg):
    """Check if a file is a FASTQ file based on its first four lines."""
    if input_arg.endswith('.fastq') or input_arg.endswith('.fq'):
        return True
    else:
        return False

def is_directory_of_fastas(input_arg):
    if os.path.isdir(input_arg):
        files = os.listdir(input_arg)
        fasta_files = [f for f in files if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna')]
        if len(fasta_files)>0:
            return True
        else:
            return False
    else:
        return False

def is_datasets_zip(input_arg):
    if input_arg.endswith('.zip'):
        return True
    else:
        return False

def read_records(input_arg, input_file_type):
    seq_dict = {}
    if input_file_type == 'zip':
        with zipfile.ZipFile(input_arg, 'r') as infile:
            for f in infile.namelist():
                if is_fasta_file(f):
                    with infile.open(f) as file:
                        records = SeqIO.parse(StringIO(file.read().decode()), 'fasta')
                        for record in records:
                            if record.description not in seq_dict:
                                seq_dict[record.description] = [str(record.seq)]
                            else:
                                seq_dict[record.description].append(str(record.seq))
    
    elif input_file_type == 'directory':
        for filename in os.listdir(input_arg):
            file_path = os.path.join(input_arg, filename)
            if os.path.isfile(file_path) and is_fasta_file(filename):
                records = SeqIO.parse(file_path, 'fasta')
                for record in records:
                    if record.description not in seq_dict:
                        seq_dict[record.description] = [str(record.seq)]
                    else:
                        seq_dict[record.description].append(str(record.seq))

    if input_file_type == 'fasta':
        for record in SeqIO.parse(input_arg, "fasta"):
            # Use the record.id as the key and the sequence as the value
            seq_dict[record.description] = [str(record.seq)]
            
    if input_file_type == 'fastq':
        for record in SeqIO.parse(input_arg, "fastq"):
            # Use the record.id as the key and the sequence as the value
            seq_dict[record.description] = [str(record.seq)]
            
    return seq_dict

def n_filter_type(value):
    # Check if the value is "any"
    if value == "any":
        return value
    try:
        # First, try to convert the value to a float
        float_value = float(value)
        
        # Check if it's actually an integer (no decimal part)
        if float_value.is_integer():
            return int(float_value)  # Return as int if it's an integer
        else:
            return float_value  # Return as float if it's a true float
    except ValueError:
        # If conversion fails, the input is neither a float nor an integer nor "any"
        raise argparse.ArgumentTypeError(f"Invalid n_filter value: {value}")

def contains_degenerate_characters(sequence):
    allowed_chars = {'A', 'T', 'C', 'G'}
    for char in sequence.upper():
        if char not in allowed_chars:
            return True  # Found a degenerate character
    return False  # No degenerate characters found

def calculate_gc_content(sequence):
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    if len(sequence) == 0:  # Avoid division by zero for empty sequences
        return 0
    gc_content = (gc_count / len(sequence))
    return gc_content

def header_contains_in_strings(header, filter_strings):
    # Convert header to lower case for case-insensitive comparison
    header_lower = header.lower()
    flag=True
    for filter_string in filter_strings:
        # Convert filter string to lower case too
        filter_string_lower = filter_string.lower()
        if filter_string_lower not in header_lower:
            flag = False  # Found a filter string not in header
    return flag  # All filter strings found in header

def header_contains_out_strings(header, filter_strings):
    # Convert header to lower case for case-insensitive comparison
    header_lower = header.lower()
    flag=False
    for filter_string in filter_strings:
        # Convert filter string to lower case too
        filter_string_lower = filter_string.lower()
        if filter_string_lower in header_lower:
            flag = True  # Found a filter string not in header
    return flag  # All filter strings found in header


def compare_seq_dicts(original_seq_dict, modified_seq_dict):
    added_keys = set(modified_seq_dict.keys()) - set(original_seq_dict.keys())
    removed_keys = set(original_seq_dict.keys()) - set(modified_seq_dict.keys())
    changed_sequences = {key: (original_seq_dict[key], modified_seq_dict[key])
                         for key in original_seq_dict
                         if key in modified_seq_dict and original_seq_dict[key] != modified_seq_dict[key]}
    
    return added_keys, removed_keys, changed_sequences


def check_identical_values(dictionary):
    refseq_leaders = [">NM_", ">NP_", ">NR_", ">NG_", ">NT_", ">NW_", ">NZ_", ">NC_", ">XM_", ">XP_", ">XR_"]
    
    # Initialize a dictionary to hold the chosen sequences and their corresponding keys
    chosen_sequences = {}
    
    # Iterate over each key, value pair in the dictionary
    for key, sequences in dictionary.items():
        for seq in sequences:
            # Convert sequence to a hashable type for comparison and storage
            seq_tuple = tuple(seq)
            
            # Check if this sequence has already been processed
            if seq_tuple in chosen_sequences:
                existing_key = chosen_sequences[seq_tuple]
                # Determine if the current key has a higher priority based on refseq_leaders
                if any(key.startswith(prefix) for prefix in refseq_leaders) and not any(existing_key.startswith(prefix) for prefix in refseq_leaders):
                    # Replace the existing key with the current key if it has higher priority
                    chosen_sequences[seq_tuple] = key
            else:
                # Add the sequence and its key to the chosen_sequences if not already present
                chosen_sequences[seq_tuple] = key
    
    # Reconstruct the dictionary to match the original format
    result_dict = {}
    for seq_tuple, key in chosen_sequences.items():
        if key in result_dict:
            result_dict[key].append(''.join(seq_tuple))
        else:
            result_dict[key] = [''.join(seq_tuple)]
    
    return result_dict



def compare_dict_keys(dict1, dict2):
    diff1 = set(dict1.keys()) - set(dict2.keys())
    diff2 = set(dict2.keys()) - set(dict1.keys())
    return list(diff1.union(diff2))


parser = argparse.ArgumentParser("QC Check - provide datasets zip, a multi fasta, or a directory of fasta files. Perform QC as desired. Output one of the three formats")

parser.add_argument('-i', '--input', required=True, help="NCBI Datasets zip file, multifasta file, or directory of fasta files to read in")
parser.add_argument('-n', '--n_filter', required=False, default=False, help='input criteria for filtering by N characters in sequence. If a float, filters by that percent of Ns in sequence. If an integer, filters any sequences with that many consecutive Ns. If "any", filter records containing any number of N')
parser.add_argument('-d', '--degen_filter', required=False, action='store_true', help='flag if all records containing any degen characters should be filtered')
parser.add_argument('--gc_lower', required=False, default=0.0, type=float, help='lower bound of GC content to be kept in sequences')
parser.add_argument('--gc_upper', required=False, default=1.0, type=float, help='upper bound of GC content to be kept in sequences')
parser.add_argument('-bm', '--bin_method', required=False, default=False, help='input criteria for filtering by N characters in sequence. Either "median" if binning should be done by median sequence length, or specify an integer length')
parser.add_argument('-bp', '--bin_percentile',required=False, type=float, default=.1, help='float value for  percent difference in length to cut off binning (+- 10 percent by default)')
parser.add_argument('-si', '--strings_to_include', required=False, nargs='*', help='An optional space-separated list of strings to filter inclusively')
parser.add_argument('-se', '--strings_to_exclude', required=False, nargs='*', help='An optional space-separated list of strings to filter exclusively')
parser.add_argument('-on', '--output_name', required=True, type=str, help='What name should the output have')
parser.add_argument('-of', '--output_format', choices=['fasta', 'directory', 'zip'], help='Output format: fasta, directory, or zip.')
parser.add_argument('-sd', '--save_drops', required=False,default=False, action='store_true',help='true or false to write out all filtered sequences')
parser.add_argument('-cd', '--check_duplicates', required=False,default=False, action='store_true', help='drop duplicate squences with preference for refseq sequences remaining')

args = parser.parse_args()

if args.save_drops:
    drop_dict = {}

# Determine input type and read in sequence records
input_file_type = None
if is_fasta_file(args.input):
    input_file_type = "fasta"

elif is_fastq_file(args.input):
    input_file_type = "fastq"

elif is_directory_of_fastas(args.input): # this form will be filename by list of seqs
    input_file_type = "directory"
    
elif is_datasets_zip(args.input): # this form will be filename by list of seqs
    input_file_type = "zip"
    
else:
    raise Exception("An error occurred: Invalid input")
    
seq_dict = read_records(args.input, input_file_type)
original_seq_dict = copy.deepcopy(seq_dict)

if args.check_duplicates:
    print("Filtering duplicate sequences")
    filtered_dict = check_identical_values(seq_dict)
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'Duplicate_record'
    seq_dict = filtered_dict

if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")
    
# Apply QC steps
n_filter = args.n_filter
# N content
if n_filter:
    print("Filtering N character sequences")
    n_content_value = n_filter_type(n_filter)
    filtered_dict = {}
    if isinstance(n_content_value, int):
        for header, sequences in seq_dict.items():
            # Filter each sequence in the list based on the criterion
            filtered_sequences = [seq for seq in sequences if 'N' * n_content_value not in seq.upper()]
            # Only add to the filtered dictionary if there are remaining sequences
            if filtered_sequences:
                filtered_dict[header] = filtered_sequences

            
    elif isinstance(n_content_value, float):
        for header, sequences in seq_dict.items():
            # Filter each sequence in the list based on the criterion
            filtered_sequences = [seq for seq in sequences if seq.upper().count('N') / len(seq) < n_content_value]
            # Only add to the filtered dictionary if there are remaining sequences
            if filtered_sequences:
                filtered_dict[header] = filtered_sequences

        
    elif n_content_value == 'any':

        for header, sequences in seq_dict.items():
            # Filter each sequence in the list based on the criterion
            filtered_sequences = [seq for seq in sequences if 'N' not in seq.upper()]
            # Only add to the filtered dictionary if there are remaining sequences
            if filtered_sequences:
                filtered_dict[header] = filtered_sequences

    else:
        raise Exception("An error occurred: Invalid specification for --n_filter")
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'N_filter'
    seq_dict = filtered_dict


if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")
# Check for any degen characters and report records with degen characters

if args.degen_filter:
    filtered_dict = {}
    print('Filtering degenerate characters')
    for header, sequences in seq_dict.items():
        # Filter each sequence in the list based on the criterion
        filtered_sequences = [seq for seq in sequences if not contains_degenerate_characters(seq)]
        # Only add to the filtered dictionary if there are remaining sequences
        if filtered_sequences:
            filtered_dict[header] = filtered_sequences
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'Degenerate_char_filter'
    seq_dict = filtered_dict       


if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")
# Filter by keyword list

if args.gc_lower > 0.0 or args.gc_upper < 1.0:
    filtered_dict = {}
    print('Filtering GC characters')
    for header, sequences in seq_dict.items():
        # Filter each sequence in the list based on the criterion
        filtered_sequences = [seq for seq in sequences if args.gc_lower <= calculate_gc_content(seq) <= args.gc_upper]
        # Only add to the filtered dictionary if there are remaining sequences
        if filtered_sequences:
            filtered_dict[header] = filtered_sequences
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'GC_filter'
    seq_dict = filtered_dict    

if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")

if args.strings_to_include:
    print('Filtering include strings')

    filtered_dict = {header: seq_dict[header] for header in seq_dict if header_contains_in_strings(header, args.strings_to_include)}

    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))   
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'String_inclusion_filter'         
    seq_dict = filtered_dict 


if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")

if args.strings_to_exclude:
    print('Filtering exclude strings')
    
    filtered_dict = {header: seq_dict[header] for header in seq_dict if not header_contains_out_strings(header, args.strings_to_exclude)}
                  
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'String_exclusion_filter'
    seq_dict = filtered_dict

if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")
    
# Filter by binning sequence len
if args.bin_method:
    print('Filtering by binning seq lengths')
    input_type = determine_type(args.bin_method)
    filtered_dict = {}
    if input_type == 'string':
        if args.bin_method == 'median':
            sequence_lengths = [len(sequence[0]) for header, sequence in seq_dict.items()]
            sorted_lengths = sorted(sequence_lengths)
            num_sequences = len(sorted_lengths)
            if num_sequences % 2 == 0:
                # If even number of sequences, average the middle two lengths
                median_length = (sorted_lengths[num_sequences // 2 - 1] + sorted_lengths[num_sequences // 2]) / 2
            else:
                # If odd number of sequences, take the middle length
                median_length = sorted_lengths[num_sequences // 2]
            
            upper_bound=(1+args.bin_percentile)*median_length
            lower_bound=(1-args.bin_percentile)*median_length
            print("Median: "+str(median_length))
            print("Upper bound: "+str(upper_bound))
            print("Lower bound: "+str(lower_bound))
            for header, sequence in seq_dict.items():
                if lower_bound<len(sequence[0])<upper_bound:
                    filtered_dict[header] = sequence
                            
        else:
            raise Exception("An error occurred: Invalid string specification for --bin_method")
    
    elif input_type == 'int':
        upper_bound=(1+args.bin_percentile)*int(args.bin_method)
        lower_bound=(1-args.bin_percentile)*int(args.bin_method)
        print("Median: "+str(int(args.bin_method)))
        print("Upper bound: "+str(upper_bound))
        print("Lower bound: "+str(lower_bound))
        for header, sequence in seq_dict.items():
            if lower_bound<len(sequence[0])<upper_bound:
                filtered_dict[header] = sequence
                                
    else:
        raise Exception("An error occurred: Invalid specification for --bin_method")
    print('Before filtering: '+str(len(seq_dict)))
    print('After filtering: '+str(len(filtered_dict)))
    if args.save_drops:
        dropped_keys = compare_dict_keys(seq_dict, filtered_dict)
        for key in dropped_keys:
            drop_dict[key.split(' ')[0]] = 'Binning_filter'
    seq_dict = filtered_dict 


if len(seq_dict) == 0:
    raise Exception("An error occurred: Filtered everything out!")

if args.output_format == 'fasta':
    seq_records = []
    for header, sequences in seq_dict.items():
        for sequence in sequences:
            # Create a SeqRecord object for each sequence
            seq_record = SeqRecord(Seq(sequence),
                                   id=header, 
                                   description='')  
            seq_records.append(seq_record)
    
    # Write all SeqRecord objects to a FASTA file
    with open(args.output_name, 'w') as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")

if args.output_format == 'directory':
    os.makedirs(args.output_name, exist_ok=True)
    
    for header, sequences in seq_dict.items():
        # Extract the accession number as the first part of the header
        accession = header.split(' ')[0]
        # Define the output path for each FASTA file
        output_path = os.path.join(args.output_name, f"{accession}.fasta")
        
        seq_records = []
        for sequence in sequences:
            # Create a SeqRecord object for each sequence
            seq_record = SeqRecord(Seq(sequence),
                                   id=header, 
                                   description='')  
            seq_records.append(seq_record)
        
        # Write the SeqRecord object to a FASTA file
        with open(output_path, 'w') as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")

if args.output_format == 'zip':
    base_dir = args.output_name
    ncbi_dataset_dir = os.path.join(base_dir, 'ncbi_dataset')
    data_dir = os.path.join(ncbi_dataset_dir, 'data')
    
    # Ensure base directory is clean
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
    os.makedirs(data_dir, exist_ok=True)
    
    for accession, sequences in seq_dict.items():
        accession_dir = os.path.join(data_dir, accession.split(' ')[0])
        os.makedirs(accession_dir, exist_ok=True)
        
        # Assuming each value in seq_dict is a list of sequences for the accession
        for i, sequence in enumerate(sequences):
            # If there are multiple sequences per accession, enumerate them
            fasta_filename = f"{accession.split(' ')[0]}_{i}.fna" if len(sequences) > 1 else f"{accession.split(' ')[0]}.fna"
            fasta_file_path = os.path.join(accession_dir, fasta_filename)
            
            with open(fasta_file_path, 'w') as fasta_file:
                fasta_file.write(f">{accession}\n{sequence}\n")
    
    # Zip the directory
    zip_file_path = f"{args.output_name}.zip"
    with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, start=base_dir)
                zipf.write(file_path, arcname)
    
    # Remove the temporary directory structure
    shutil.rmtree(base_dir)    


if args.save_drops:
    _, removed_keys, changed_sequences = compare_seq_dicts(original_seq_dict, seq_dict)
    seq_records = []
    for key in removed_keys:
        sequence = original_seq_dict[key]
        seq_record = SeqRecord(Seq(sequence[0]), id=key, description="")
        seq_records.append(seq_record)

    # Write to FASTA
    with open("REMOVED_RECORDS_"+args.output_name+'.fasta', 'w') as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")
    
    with open("REMOVED_RECORDS_"+args.output_name+'_meta.tsv', 'w') as file:
        # Write the header
        file.write("accession\tfilter_criteria\n")
        # Write the dictionary content
        for key, value in drop_dict.items():
            file.write(f"{key}\t{value}\n")





















