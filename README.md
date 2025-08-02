# CURE2025_Comparative_massSpec

## Overview

This repository contains a comprehensive Jupyter notebook for **Processing Mass Spectrometry Proteomics Data**. The notebook provides an in-depth tutorial on mass spectrometry fundamentals and hands-on analysis of proteomics datasets from two different organisms:

- **Soybean (*Glycine max*)** - Label-free quantification (LFQ) analysis of phosphate stress response
- **Human (*Homo sapiens*)** - SILAC-labeled analysis of hypoxia-induced pulmonary fibrosis

### What You'll Learn

- Fundamentals of mass spectrometry and proteomics workflows
- How to process raw mass spectrometry data using SAGE (Spectral Alignment Guided Engine)
- Peptide identification and protein quantification techniques
- Statistical analysis of differential protein expression
- Data visualization and interpretation of proteomics results

## Prerequisites

- Basic knowledge of Python programming
- Understanding of biological concepts (proteins, enzymes, cellular processes)
- Familiarity with command-line operations (helpful but not required)

## Environment Setup

### Step 1: Install Conda/Miniconda

If you don't have conda installed, download and install Miniconda from:
https://docs.conda.io/en/latest/miniconda.html

### Step 2: Create the Conda Environment

Navigate to the repository directory and create the environment using the provided environment file:

```bash
# Navigate to the project directory
cd CURE2025_Comparative_massSpec

# Create the conda environment
conda env create -f src/data/environment.yml

# Activate the environment
conda activate mass_spec_proteomics
```

### Step 3: Install SAGE (Optional)

The notebook uses SAGE for proteomics analysis. If you want to run the SAGE portions:

1. Download SAGE from: https://github.com/lazear/sage
2. Follow the installation instructions for your operating system
3. Ensure SAGE is available in your system PATH

### Step 4: Launch Jupyter Notebook

If you are on CyVerse:
1. Open the file browser denoted by the folder icon on the upper left corner. 
2. Navigate to the following directory `src/data/` and then double click the jupyter notebook named `Processing_Mass_Spectrometry_Proteomics_Data.ipynb`.  

If you are on a command line:
```bash
# Make sure the environment is activated
conda activate mass_spec_proteomics

# Start Jupyter notebook
jupyter notebook
```

This will open your default web browser with the Jupyter interface.

## Jupyter Notebook Tutorial

### Getting Started with Jupyter

1. **Opening the Notebook**
   - Navigate to `src/data/Processing_Mass_Spectrometry_Proteomics_Data.ipynb`
   - Click on the file to open it

2. **Understanding Cell Types**
   - **Markdown cells**: Contain text, explanations, and formatted content
   - **Code cells**: Contain Python code that can be executed

### Navigation Basics

- **Moving between cells**: Use arrow keys or click on cells
- **Cell selection**: Click on a cell to select it (blue border = selected, green border = edit mode)
- **Enter edit mode**: Press `Enter` or double-click on a cell
- **Exit edit mode**: Press `Esc`

### Running Cells

- **Run current cell**: Press `Shift + Enter` (runs cell and moves to next)
- **Run current cell (stay)**: Press `Ctrl + Enter` (runs cell but stays on current cell)
- **Run all cells**: Menu > Cell > Run All
- **Run cells above**: Menu > Cell > Run All Above
- **Run cells below**: Menu > Cell > Run All Below

### Keyboard Shortcuts (Command Mode)

When a cell is selected but not in edit mode (blue border):

- `A` - Insert cell above
- `B` - Insert cell below
- `DD` - Delete cell (press D twice)
- `M` - Change cell to Markdown
- `Y` - Change cell to Code
- `Z` - Undo cell deletion
- `Shift + M` - Merge selected cells

### Working with Code Cells

1. **Executing Code**
   ```python
   # Example: Run this cell to see output
   print("Hello, proteomics!")
   ```

2. **Variables persist**: Variables created in one cell are available in later cells

3. **Viewing outputs**: Results appear below the cell after execution

4. **Handling errors**: If a cell produces an error, read the traceback to understand the issue

### Best Practices for This Notebook

1. **Run cells in order**: The notebook is designed to be run sequentially from top to bottom

2. **Read the markdown cells**: They contain important explanations and context

3. **Don't skip imports**: Early cells contain necessary library imports

4. **Check file paths**: Ensure data files are in the correct locations

5. **Monitor memory usage**: Large datasets may require significant RAM

### Troubleshooting Common Issues

1. **Kernel issues**
   - If the kernel becomes unresponsive: Kernel > Restart
   - To restart and run all cells: Kernel > Restart & Run All

2. **Import errors**
   - Ensure the conda environment is activated
   - Check that all required packages are installed

3. **File not found errors**
   - Verify file paths in the notebook match your directory structure
   - Check that data files have been downloaded correctly

4. **Memory errors**
   - Close other applications to free up RAM
   - Consider running on a machine with more memory

### Saving Your Work

- **Auto-save**: Jupyter automatically saves periodically
- **Manual save**: Press `Ctrl + S` or use File > Save and Checkpoint
- **Export options**: File > Download as > (HTML, PDF, etc.)

### Advanced Features

1. **Magic commands**
   ```python
   %matplotlib inline  # Display plots inline
   %%time             # Time execution of a cell
   %pwd               # Show current directory
   ```

2. **Getting help**
   ```python
   help(function_name)  # Get help for a function
   function_name?       # Quick help in Jupyter
   ```

3. **Tab completion**: Press `Tab` while typing to see available options

## Notebook Structure

The notebook is organized into sections:

1. **Section I**: Mass Spectrometry Fundamentals
   - Theory and principles
   - Top-down vs bottom-up proteomics
   - Instrumentation overview

2. **Section II**: Hands-on Analysis
   - Data collection and quality control
   - SAGE configuration and execution
   - Statistical analysis and visualization

## Data Requirements

The notebook expects data in specific locations. Ensure you have:
- Raw mass spectrometry files (.mzML format)
- Protein databases (FASTA format)
- Appropriate directory structure as shown in the notebook

## Support

If you encounter issues:
1. Check the troubleshooting section above
2. Verify your environment setup
3. Ensure all dependencies are correctly installed
4. Review error messages carefully for clues

For additional help with Jupyter notebooks, visit: https://jupyter-notebook.readthedocs.io/



