## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/sgrolab/cellsimilaritymodel.git
cd cellsimilaritymodel
```

### 2. Install uv
If you don't already have `uv` installed, install it using one of the following methods:

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows:**
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (using pip):**
```bash
pip install uv
```

### 3. Create Virtual Environment and Install Dependencies
```bash
# Create a virtual environment with uv
uv venv

# Activate the virtual environment
# On macOS/Linux:
source .venv/bin/activate
# On Windows:
.venv\Scripts\activate

# Install project dependencies
uv pip install -r requirements.txt
```

### 4. Configure Project Settings
Create a configuration file to specify your data directory:
```bash
# Create the utils directory if it doesn't exist
mkdir -p utils

# Create the config.py file
cat > utils/config.py << EOF
"""
Configuration file for project paths and settings
"""
import os
from pathlib import Path

# Set the project directory as the home directory for data
PROJECT_DIR = Path("/path/to/your/data/directory")

# Create data directory if it doesn't exist
PROJECT_DIR.mkdir(parents=True, exist_ok=True)

EOF
```

**Important:** Edit `utils/config.py` and replace `/path/to/your/data/directory` with your actual data directory path.

### 5. Verify Installation
```bash
# Test that the environment is set up correctly
python -c "from utils.config import PROJECT_DIR; print(f'Project directory: {PROJECT_DIR}')"
```

You should see your configured project directory path printed to the console.