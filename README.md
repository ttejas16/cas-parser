# Project Setup Instructions

## Prerequisites
- Python 3.6 or higher installed on your system
- pip package manager

## Setting up Virtual Environment

### 1. Create Virtual Environment
```bash
python -m venv venv
```

### 2. Activate Virtual Environment

**On Windows:**
```bash
venv\Scripts\activate
```

**On macOS/Linux:**
```bash
source venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Verify Installation
```bash
pip list
```

## Running the Project
Once the virtual environment is activated and dependencies are installed, you can run the project:
```bash
python main.py
```

## Deactivating Virtual Environment
When you're done working on the project:
```bash
deactivate
```

## Troubleshooting

### If you encounter permission errors:
- On Windows: Run command prompt as administrator
- On macOS/Linux: Use `sudo` if necessary, but prefer using `--user` flag

### If virtual environment creation fails:
```bash
# Try using python3 explicitly
python3 -m venv venv

# Or if you have multiple Python versions
py -3 -m venv venv
```

### If pip install fails:
```bash
# Upgrade pip first
python -m pip install --upgrade pip

# Then install requirements
pip install -r requirements.txt
```

## Notes
- Always activate the virtual environment before working on the project
- The `venv` folder should not be committed to version control (add to `.gitignore`)
- If you add new packages, update requirements.txt with: `pip freeze > requirements.txt`