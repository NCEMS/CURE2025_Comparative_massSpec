import glob
import multiprocessing
import requests
import os, sys
import subprocess

##############################################################################
def download_ftp_file_wget(ftp_url, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = ["wget", "-c", ftp_url, "-P", output_dir]
    subprocess.run(cmd, check=True)
##############################################################################

##############################################################################
def convert_RAW(file_path: str):
    """
    Convert a .raw file to a .mgf and .mzML format using ThermoRawFileParser.
    """
    try:
        print(f"Converting {file_path} --> .mgf")
        if not os.path.exists(file_path.replace('.raw', '.mgf')):
            cmd = ['ThermoRawFileParser', '-i', file_path, '-f', '0']
            subprocess.run(cmd, check=True)
            print(f"Finished converting {file_path} --> .mgf")
        else:
            print(f"MGF file already exists: {file_path.replace('.raw', '.mgf')}")
    except Exception as e:
        print(f"Error converting {file_path} to .mgf: {e}")


    try:
        print(f"Converting {file_path} --> .mzML")
        if not os.path.exists(file_path.replace('.raw', '.mzML')):
            cmd = ['ThermoRawFileParser', '-i', file_path, '-f', '1']
            subprocess.run(cmd, check=True)
            print(f"Finished converting {file_path} --> .mzML")
        else:
            print(f"MZML file already exists: {file_path.replace('.raw', '.mzML')}")
    except Exception as e:
        print(f"Error converting {file_path} to .mzML: {e}")
##############################################################################

##############################################################################
def download_pride_files(pxd_id, file_type="", output_dir="raw_data"):
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd_id}/files"
    response = requests.get(url)
    response.raise_for_status()

    file_list = response.json()
    print(f'Number of files: {len(file_list)}')

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Filter and download matching files
    file_list = [(f, output_dir) for f in file_list if f["fileName"].endswith(file_type)]

    # Use multiprocessing to process files in parallel
    with multiprocessing.Pool() as pool:
        pool.map(process_file, file_list)
    ####################################################

    ####################################################
    # Launch the RunAssessor command
    print(f'\nLaunching Run Assessor')
    runAssessor_output_dir = os.path.join(output_dir, 'runAssessor')
    runAssessor_outfile = os.path.join(runAssessor_output_dir, 'study_metadata.json')
    if os.path.exists(runAssessor_outfile):
        print(f"RunAssessor output already exists: {runAssessor_outfile}")
    else:
        cmd = f"python src/data/mzML_assessor.py --inpath {output_dir} --outpath {runAssessor_output_dir}"
        print(f"{'#'*50}\nRunning command: {cmd}")
        cmd_split = cmd.split()
        print(cmd_split)
        result = subprocess.run(cmd_split, capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.returncode != 0:
            raise SyntaxError(f'Error occurred while running mzML_assessor:\n{result.stderr}')
    ####################################################

    print("Downloads complete.")
##############################################################################

##############################################################################
def process_file(args):
    f, output_dir = args
    # Check if the file type matches the specified type
    file_name = f["fileName"]
    print(f"File name: {file_name}")

    # Check if the file already exists
    outfile = os.path.join(output_dir, file_name)
    print(f"Output file path: {outfile}")
    if os.path.exists(outfile):
        print(f"File {outfile} already exists, skipping download.")
    else:
        download_url = f['publicFileLocations']
        FTP_url = ''
        if isinstance(download_url, list):
            for url in download_url:
                if "FTP" in url['name']:
                    FTP_url = url['value']
                    break
        print(f"FTP URL: {FTP_url}")
        if isinstance(FTP_url, list):
            raise ValueError("Download URL is a list, expected a single URL. No FTP protocol is avaiable for this file.")
        print(f"Downloading {file_name}...")
        download_ftp_file_wget(FTP_url, output_dir)
        print(f"Downloaded {file_name} to {output_dir}")

    # Convert the downloaded .raw file to .mgf and .mzML formats
    convert_RAW(outfile)
    print(f"Converted {file_name} to .mgf and .mzML formats.")
##############################################################################

##############################################################################
import argparse
import pandas as pd

# Example usage:
parser = argparse.ArgumentParser(description="Download PRIDE files")
parser.add_argument("--input", help="Either a single PXD ID string or a .csv file where the column named PXD is used")
parser.add_argument("--output_dir", help="Output directory")
args = parser.parse_args()

PXD = args.input
outdir = args.output_dir

## check if the PXD argument is a single ID or a CSV file
if PXD.endswith(".csv"):
    print(f"Detected CSV file: {PXD}")
    df = pd.read_csv(PXD)
    if 'PXD' not in df.columns:
        raise ValueError("CSV file must contain a 'PXD' column.")
    PXDs = df['PXD'].unique().tolist()
else:
    print(f"Detected single PXD ID: {PXD}")
    PXDs = [PXD]

print(f'Downloading PXDs {PXDs} {len(PXDs)}')
print(f'Output directory: {outdir}')

for PXD in PXDs:
    print(f"Processing PXD: {PXD}")
    PXD_outpath = os.path.join(outdir, PXD)

    ## Download files of type .raw
    download_pride_files(PXD, file_type=".raw", output_dir=PXD_outpath)

print('NORMAL TERMINATION')
##############################################################################