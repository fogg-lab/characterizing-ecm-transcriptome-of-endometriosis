# Dataset pathing
## Developer configuration file
Make sure you have a file, "dev_paths.txt", in the root project directory laid out as follows:
```
root/path/for/storing/data (e.g. D:/unified_TCGA_GTEx)
<newline>
```

## Directory verification script
Run
```
python verify_dirs.py
```
from the project root. This will verify that all needed files/directories exist.
