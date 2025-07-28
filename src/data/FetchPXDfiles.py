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
def download_pride_files(pxd_id, file_type="", output_dir="raw_data"):
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd_id}/files"
    response = requests.get(url)
    response.raise_for_status()

    file_list = response.json()
    print(f'Number of files: {len(file_list)}')

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Filter and download matching files
    for f in file_list:
        #print(f)
      
        # Check if the file type matches the specified type
        file_name = f["fileName"]
        print(f"File name: {file_name}")
        if file_type != "":
            if not file_name.endswith(file_type):
                print(f"Skipping {file_name} (not a {file_type} file)")
                continue

        download_url = f['publicFileLocations']
        FTP_url = ''
        if isinstance(download_url, list):
            for url in download_url:
                #print(url)
                if "FTP" in url['name']:
                    FTP_url = url['value']
                    break
        print(f"FTP URL: {FTP_url}")

        # Check if the file name ends with the specified file type
        if isinstance(FTP_url, list):
            raise ValueError("Download URL is a list, expected a single URL. No FTP protocol is avaiable for this file.")

        
        print(f"Downloading {file_name}...")
        download_ftp_file_wget(FTP_url, output_dir)
        print(f"Downloaded {file_name} to {output_dir}")
        
        # Stop after downloading the first file
        break

    print("Downloads complete.")
##############################################################################

# Example usage:
PXD = sys.argv[1]
outdir = sys.argv[2]
print(f'Downloading PXD {PXD}')
print(f'Output directory: {outdir}')
# download_pride_files(PXD, file_type=".raw", output_dir=outdir)
download_pride_files(PXD, file_type="", output_dir=outdir)

