# README — Open and Run a Jupyter Notebook in VS Code (Mac & Windows/WSL)

This guide walks students through installing **Visual Studio Code**, the **Python** + **Jupyter** extensions, setting up **Python**, and launching your `.ipynb` notebook. Windows users will install and use **WSL (Windows Subsystem for Linux)** for a smoother, Linux-like experience.

---

## TL;DR (what students do)

1. **Install VS Code** and the **Python** and **Jupyter** extensions. ([Visual Studio Code][1])
2. **Install Python** (macOS: from python.org; Windows: via WSL). ([Python.org][2])
3. **Create a virtual environment** and **install Jupyter** (`pip install jupyter ipykernel`).
4. **Open the course folder in VS Code** → open the `.ipynb` → **Select Kernel** (pick your venv) → **Run**.

---

## What you’ll install

* **Visual Studio Code** (editor) and extensions: **Python** + **Jupyter**. ([Visual Studio Code][1])
* **Python 3.10+** (3.11+ recommended). ([Python.org][3])
* **WSL** (Windows only) + **Ubuntu** distro, plus Python tools in WSL. ([Microsoft Learn][4])

---

## macOS setup

### 1) Install VS Code

* Download and install VS Code for macOS. Launch it once. ([Visual Studio Code][1])

### 2) Install the extensions

Open VS Code → **Extensions** (left sidebar) → search and install:

* **Python** (Microsoft)
* **Jupyter** (Microsoft)
  These provide Notebook support and Python tooling. ([Visual Studio Code][5])

### 3) Install Python

* Easiest: download the macOS installer from **python.org** and run it (ensure “Add to PATH” is enabled if prompted). After install, `python3 --version` in **Terminal** should show a recent 3.x. ([Python.org][2])

### 4) Make a project folder and virtual environment

Open **Terminal** (Spotlight → “Terminal”) and run:

```bash
# replace "my-course" with your folder name
mkdir -p ~/my-course && cd ~/my-course

# create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# upgrade pip and install Jupyter + common basics
python -m pip install -U pip
pip install jupyter ipykernel
```

> Optional (nice to have):
> `python -m ipykernel install --user --name my-course --display-name "Python (my-course)"`
> This gives your kernel a friendly name in VS Code.

### 5) Open the notebook and run it

* In VS Code: **File → Open Folder…** → select your project folder.
* Open the `.ipynb`.
* If prompted, click **Trust**.
* Top-right of the notebook, click **Select Kernel** → choose your venv (or “Python (my-course)”).
* Press **Run** ▷ on cells.

**Troubleshooting (macOS)**

* **No kernels found**: make sure your venv is activated and `pip install jupyter ipykernel` completed without errors.
* **Wrong interpreter**: VS Code Status Bar → click the Python version → choose the interpreter from `.venv`. ([Visual Studio Code][5])

---

## Windows setup (with WSL — recommended)

We’ll run Python in **Ubuntu on WSL**, and use VS Code’s **Remote – WSL** workflow.

### 1) Install VS Code (Windows)

Download and install VS Code (User Installer). Launch it once. ([Visual Studio Code][1])

### 2) Install WSL + Ubuntu

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

* Restart when prompted.
* On first launch of “Ubuntu”, create a **username** and **password**.
* Update packages:

```bash
sudo apt update && sudo apt -y upgrade
```

WSL defaults to **WSL 2** on current Windows; you can verify with `wsl -l -v`. ([Microsoft Learn][4])

### 3) Install Python tools inside Ubuntu (WSL)

In the Ubuntu terminal:

```bash
sudo apt install -y python3 python3-venv python3-pip
```

### 4) Install the VS Code extensions

Open VS Code (Windows) → **Extensions**:

* **Remote – WSL** (Microsoft)
* **Python** (Microsoft)
* **Jupyter** (Microsoft)
  This enables opening Linux folders and running Notebooks from WSL transparently. ([Visual Studio Code][6])

### 5) Open the project folder in WSL

Option A (from WSL terminal):

```bash
# create a project folder (or cd into your cloned repo)
mkdir -p ~/my-course && cd ~/my-course

# launch VS Code connected to WSL in this folder
code .
```

Option B (from VS Code): Click the green remote button in the lower-left (><) → **WSL: Open Folder** → choose your Linux path (e.g., `/home/<you>/my-course`). ([Visual Studio Code][7])

### 6) Create a venv and install Jupyter (still inside WSL)

In VS Code’s integrated terminal (it should say “WSL: Ubuntu” in the bottom-left):

```bash
python3 -m venv .venv
source .venv/bin/activate

python -m pip install -U pip
pip install jupyter ipykernel
python -m ipykernel install --user --name my-course --display-name "Python (my-course)"
```

### 7) Open and run the notebook

* In VS Code (WSL window), open the `.ipynb`.
* **Select Kernel** → choose **Python (my-course)** (or the `.venv` interpreter).
* Run cells.

**Troubleshooting (Windows/WSL)**

* **`code .` not found**: Install the **Remote – WSL** extension and reopen your WSL terminal; the `code` command is added by the extension. ([Visual Studio Code][6])
* **Notebook won’t run**: Confirm the **active kernel** is your WSL venv and that `jupyter` is installed there. ([Visual Studio Code][8])
* **WSL version check**: `wsl -l -v`; use `wsl --set-version <Distro> 2` if needed. ([Microsoft Learn][4])

---

## (Optional) Native Windows (without WSL)

If you must run natively:

1. Install **Python for Windows** from python.org (enable “Add Python to PATH”). ([Python.org][9])
2. In **Command Prompt** (or PowerShell):

```bat
mkdir %USERPROFILE%\my-course
cd %USERPROFILE%\my-course
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install -U pip
pip install jupyter ipykernel
python -m ipykernel install --user --name my-course --display-name "Python (my-course)"
```

3. Open the folder in VS Code, open `.ipynb`, **Select Kernel** → your venv. ([Visual Studio Code][5])

> Note: WSL tends to be more consistent across student machines; prefer WSL unless you have a reason to stay native. ([Visual Studio Code][6])

---

## How to launch the notebook (any OS)

1. **Open the folder**, not just the file: **File → Open Folder…** (select the course folder).
2. Open the `.ipynb`.
3. **Select Kernel** (upper-right) → choose your venv (“Python (my-course)” or `.venv`).
4. **Run All** or run cells one by one.
5. If VS Code asks to **Trust the notebook**, click **Trust**.
6. If you get a kernel error, verify `jupyter` and `ipykernel` are installed **in that environment**, and re-select the interpreter. ([Visual Studio Code][8])

---

## Extension recap (what to search for in the Extensions view)

* **Python** — Microsoft (`ms-python.python`)
* **Jupyter** — Microsoft (`ms-toolsai.jupyter`)
* (Optional) **Jupyter Renderers** for rich outputs (`ms-toolsai.jupyter-renderers`)
* (Windows) **Remote – WSL** (`ms-vscode-remote.remote-wsl`)
  Install and manage from the Extensions view (`Ctrl+Shift+X` / `⌘⇧X`). ([Visual Studio Code][10], [Visual Studio Marketplace][11])

---

## Tips

* **Keep environments per project**: one venv per course folder avoids dependency conflicts.
* **Use the Command Palette** (`Ctrl+Shift+P` / `⌘⇧P`):

  * “**Python: Select Interpreter**” to switch venvs.
  * “**Jupyter: Select Kernel**” for notebooks. ([Visual Studio Code][5])
* **Autosave**: consider enabling **File → Auto Save** in VS Code.
* **Git**: if you’re using Git, open the **folder** as a repo and commit your work regularly.

---

## References

* **Download VS Code**; **Jupyter in VS Code**; **Python in VS Code**; **Extensions Marketplace**. ([Visual Studio Code][1])
* **Install WSL / Use WSL with VS Code**. ([Microsoft Learn][4], [Visual Studio Code][6])
* **Download Python** (macOS/Windows). ([Python.org][2])

---

### (Instructor note)

If you’re distributing a repo, add a short **`SETUP.md`** with the exact package list (or a `requirements.txt`) and remind students to **open the folder**, then **select the correct kernel** before running.

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

