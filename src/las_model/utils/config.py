import os 
from pathlib import Path 
from dotenv import load_dotenv 

load_dotenv()

PROJECT_DIR = Path(os.getenv("ROOT_DIR"))