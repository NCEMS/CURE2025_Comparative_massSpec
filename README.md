# CURE2025_Comparative_massSpec

## Overview

This repository contains a comprehensive Jupyter notebook for Processing Mass Spectrometry Proteomics Data. The notebook provides an in-depth tutorial on mass spectrometry fundamentals and hands-on analysis of proteomics datasets from two different organisms:

- **Soybean (*Glycine max*)** - Label-free quantification (LFQ) analysis of phosphate stress response
- **Human (*Homo sapiens*)** - SILAC-labeled analysis of hypoxia-induced pulmonary fibrosis

### What You'll Learn

- Fundamentals of mass spectrometry and proteomics workflows
- How to process raw mass spectrometry data using SAGE (Spectral Alignment Guided Engine)
- Peptide identification and protein quantification techniques
- Statistical analysis of differential protein expression
- Data visualization and interpretation of proteomics results

### Materials

- [MassSpec_TheoryReview.pdf](documents/MassSpec_TheoryReview.pdf): An overview of mass spectrometry in proteomics.  
- [ProjectOutline_v2.1.pdf](documents/ProjectOutline_v2.1.pdf): The project outline containing the key milestones and questions to answer. 
- [Soybean_msstats.csv](data/Soybean_msstats.csv): The processed mass spec data for the Soybean dataset. 
- [Processing_SOYBEAN_MSFragger_Mass_Spectrometry_Proteomics_Data.ipynb](src/data/Processing_SOYBEAN_MSFragger_Mass_Spectrometry_Proteomics_Data.ipynb): Example notebook to analyze the Soybean samples. 
- [Human_msstats.csv](data/Human_msstats.csv): The processed mass spec data for the Human dataset. 
- [Processing_HUMAN_MSFragger_Mass_Spectrometry_Proteomics_Data.ipynb](src/data/Processing_HUMAN_MSFragger_Mass_Spectrometry_Proteomics_Data.ipynb): Example notebook to analyze the Human samples. 

All output files from either analysis can be found in the following directories:  
[data/Soybean_processed/](data/Soybean_processed/)  
[data/Human_processed/](data/Human_processed/)  

---

## Google Collab Usage
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](
https://colab.research.google.com/github/NCEMS/CURE2025_Comparative_massSpec/blob/UndergradLevel/src/data/HUMAN_MSFragger_Mass_Spectrometry_Proteomics_Data_GoogleCollab.ipynb)
---

## Setting up your local system  

This guide walks students through installing Visual Studio Code, the Python + Jupyter extensions, setting up Python, and launching your `.ipynb` notebook. Windows users will install and use WSL (Windows Subsystem for Linux) for a smoother, Linux-like experience.

---

## macOS setup

### 1) Install VS Code

* Download and install VS Code for macOS. Launch it once. ([Visual Studio Code][1])

### 2) Install the extensions

Open VS Code → Extensions (left sidebar) → search and install:

* Python (Microsoft)
* Jupyter (Microsoft)
  These provide Notebook support and Python tooling. ([Visual Studio Code][5])

### 3) Install Miniconda

* Download and install Miniconda for macOS from [Miniconda Download](https://docs.conda.io/en/latest/miniconda.html). You can do so in the command line by the following commands:   
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
```

### 4) Clone the course repository

Open Terminal and run:

```bash
git clone -b StandAlone_v2.1.1 https://github.com/NCEMS/CURE2025_Comparative_massSpec.git
cd CURE2025_Comparative_massSpec
```

This will create a folder called `CURE2025_Comparative_massSpec` containing all course materials and notebooks.

### 5) Create the NCEMS conda environment

In the same terminal window, run:

```bash
conda create -n NCEMS python=3.11
conda activate NCEMS
conda install jupyter ipykernel pandas matplotlib seaborn scipy statsmodels
conda install -c conda-forge pyteomics
```

> Optional:
> `python -m ipykernel install --user --name NCEMS --display-name "Python (NCEMS)"`
> This gives your kernel a friendly name in VS Code.

### 6) Open the notebook and run it

* In VS Code: File → Open Folder… → select the `CURE2025_Comparative_massSpec` folder you cloned.
* Open the notebook by entering the following on the command line: 
  ```bash
  code src/data/Processing_SOYBEAN_directlfq_Mass_Spectrometry_Proteomics_Data.ipynb
  ```  
  It may install a VS Code server or start the notebook in a different window.  
* If prompted, click Trust.
* Top-right of the notebook, click Select Kernel → choose your NCEMS environment.
* Press Run ▷ on each cell in order to examine the data. 

**Troubleshooting (macOS)**

* No kernels found: make sure your venv is activated and `pip install jupyter ipykernel` completed without errors.
* Wrong interpreter: VS Code Status Bar → click the Python version → choose the interpreter from `.venv`. ([Visual Studio Code][5])

---

## Windows setup (with WSL — recommended)

We’ll run Python in Ubuntu on WSL, and use VS Code’s Remote – WSL workflow.

### 1) Install VS Code (Windows)

Download and install VS Code (User Installer). Launch it once. ([Visual Studio Code][1])

### 2) Install WSL + Ubuntu

Open PowerShell as Administrator and run:

```powershell
wsl --install -d Ubuntu
```

* Restart when prompted.
* On first launch of “Ubuntu”, create a username and password.
* Update packages:

```bash
sudo apt update && sudo apt -y upgrade
```

WSL defaults to WSL 2 on current Windows; you can verify with `wsl -l -v`. ([Microsoft Learn][4])

### 3) Install Miniconda inside Ubuntu (WSL)

* Download and install Miniconda for macOS from [Miniconda Download](https://docs.conda.io/en/latest/miniconda.html). You can do so in the command line by the following commands:   
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
```

### 4) Clone the course repository

In the Ubuntu terminal, run:

```bash
git clone -b StandAlone_v2.1.1 https://github.com/NCEMS/CURE2025_Comparative_massSpec.git
cd CURE2025_Comparative_massSpec
```

### 5) Install the VS Code extensions

Open VS Code (Windows) → Extensions:

* Remote – WSL (Microsoft)
* Python (Microsoft)
* Jupyter (Microsoft)
  This enables opening Linux folders and running Notebooks from WSL transparently. ([Visual Studio Code][6])

### 6) Open the project folder in WSL

There are two main ways to open your project folder in WSL, depending on your workflow:

**Option A: From the WSL terminal**

1. Open your Ubuntu (WSL) terminal.
2. Navigate to your project folder (if you haven't already):
  ```bash
  cd ~/CURE2025_Comparative_massSpec
  ```
3. Launch VS Code connected to WSL in this folder:
  ```bash
  code .
  ```
  This will open VS Code in Windows, but the workspace will be your Linux (WSL) environment. You should see "WSL: Ubuntu" in the bottom-left corner of VS Code.

  > **Tip:** If you get a "command not found" error for `code`, make sure you have installed the Remote – WSL extension in VS Code and restarted your WSL terminal. ([Visual Studio Code][6])

**Option B: From VS Code (Windows)**

1. Open VS Code on Windows.
2. Click the green remote button in the lower-left corner (><).
3. Select "WSL: Open Folder".
4. Browse to your Linux path, e.g., `/home/<your-username>/CURE2025_Comparative_massSpec`.
5. Click "OK" to open the folder in a WSL session.

  > **Note:** You may be prompted to install the VS Code server in WSL the first time you connect. Follow the prompts to complete setup.

**Troubleshooting:**

- If you do not see "WSL: Ubuntu" in the bottom-left, you may not be connected to WSL. Try reopening VS Code using `code .` from your WSL terminal.
- If the folder does not appear, double-check your path and permissions in WSL.
- If you cannot open the folder, ensure the Remote – WSL extension is installed and enabled in VS Code.

Once the folder is open in VS Code (WSL), you can use the integrated terminal, run notebooks, and install packages directly in your Linux environment.

### 7) Create the NCEMS conda environment and install Jupyter (inside WSL)

In VS Code’s integrated terminal (it should say “WSL: Ubuntu” in the bottom-left):

```bash
conda create -n NCEMS python=3.11
conda activate NCEMS
conda install jupyter ipykernel pandas matplotlib seaborn scipy statsmodels
conda install -c conda-forge pyteomics
```

### 8) Open and run the notebook

* Open the notebook by entering the following on the command line: 
```bash
code src/data/Processing_SOYBEAN_directlfq_Mass_Spectrometry_Proteomics_Data.ipynb
```  
It may install a VS Code server or start the notebook in a different window.  
* If prompted, click Trust.
* Top-right of the notebook, click Select Kernel → choose your NCEMS environment.
* Press Run ▷ on each cell in order to examine the data. 

**Troubleshooting (Windows/WSL)**

* `code .` not found: Install the Remote – WSL extension and reopen your WSL terminal; the `code` command is added by the extension. ([Visual Studio Code][6])
* Notebook won’t run: Confirm the active kernel is your WSL venv and that `jupyter` is installed there. ([Visual Studio Code][8])
* WSL version check: `wsl -l -v`; use `wsl --set-version <Distro> 2` if needed. ([Microsoft Learn][4])

---

## Tips

* Use a single conda environment per project: create and use the `NCEMS` environment for this course to avoid dependency conflicts.
* Use the Command Palette (`Ctrl+Shift+P` / `⌘⇧P`):

  * “Python: Select Interpreter” to switch conda envs.
  * “Jupyter: Select Kernel” for notebooks. ([Visual Studio Code][5])
* Autosave: consider enabling File → Auto Save in VS Code.
* Git: if you’re using Git, open the folder as a repo and commit your work regularly.

---

## References

* Download VS Code; Jupyter in VS Code; Python in VS Code; Extensions Marketplace. ([Visual Studio Code][1])
* Install WSL / Use WSL with VS Code. ([Microsoft Learn][4], [Visual Studio Code][6])
* Download Python (macOS/Windows). ([Python.org][2])

---

[1]: https://code.visualstudio.com/download?utm_source=chatgpt.com "Download Visual Studio Code - Mac, Linux, Windows"
[2]: https://www.python.org/downloads/?utm_source=chatgpt.com "Download Python"
[3]: https://www.python.org/?utm_source=chatgpt.com "Welcome to Python.org"
[4]: https://learn.microsoft.com/en-us/windows/wsl/install?utm_source=chatgpt.com "How to install Linux on Windows with WSL"
[5]: https://code.visualstudio.com/docs/languages/python?utm_source=chatgpt.com "Python in Visual Studio Code"
[6]: https://code.visualstudio.com/docs/remote/wsl?utm_source=chatgpt.com "Developing in WSL"
[7]: https://code.visualstudio.com/docs/remote/wsl-tutorial?utm_source=chatgpt.com "Remote development in WSL"
[8]: https://code.visualstudio.com/docs/datascience/jupyter-notebooks?utm_source=chatgpt.com "Jupyter Notebooks in VS Code"
[9]: https://www.python.org/getit/windows/?utm_source=chatgpt.com "Python Releases for Windows"
[10]: https://code.visualstudio.com/docs/configure/extensions/extension-marketplace?utm_source=chatgpt.com "Extension Marketplace"
[11]: https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter&utm_source=chatgpt.com "Jupyter Extension for Visual Studio Code"

