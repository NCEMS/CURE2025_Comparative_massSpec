# CURE2025_Comparative_massSpec

## Fetching PXD data files
Here's a complete documentation page for your **PRIDE dataset fetch script**, suitable for a README or internal workflow guide.

---

# 🧬 PRIDE PXD File Fetcher

This Python script allows you to automatically download mass spectrometry data files from the [PRIDE Archive](https://www.ebi.ac.uk/pride/) using a **PXD accession ID**. It retrieves file metadata from PRIDE's REST API and downloads selected files via **FTP using `wget`**.

---

## 📦 Features

* Automatically queries PRIDE Archive for files associated with a PXD ID
* Supports optional filtering by file type (e.g. `.mgf`, `.raw`, `.mzML`)
* Downloads files using FTP protocol with resume support (`wget -c`)
* Organizes output into a specified local directory

---

## 🚀 How to Use

### 1. **Basic Command**

```bash
python fetch_pride_files.py PXD003037 output_directory
```

This will:

* Fetch all files associated with `PXD003037`
* Download them to `output_directory`
* Use FTP protocol if available

---

## 🛠️ Dependencies

* Python ≥ 3.6
* Requires the following Python modules:

  * `requests`
  * `os`, `subprocess`, `sys` (standard library)
* External:

  * [`wget`](https://www.gnu.org/software/wget/) must be installed and accessible via `$PATH`

Install Python dependencies with:

```bash
pip install requests
```

Install `wget` on Ubuntu/Debian:

```bash
sudo apt install wget
```
---
## 🧠 Notes and Caveats

* Some files may not have FTP download locations (e.g., restricted or private datasets).
* The script uses `wget` for robustness (supports resume and large files), but could be adapted to use `urllib`, `requests`, or `aria2c`.
* Ensure your network allows passive FTP.
---

