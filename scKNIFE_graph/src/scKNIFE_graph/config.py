from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SRC_DIR.parent.parent

RAW_DIR = PROJECT_ROOT / "data" / "raw"
SEP_DIR = PROJECT_ROOT / "data" / "separated"
PRO_DIR = PROJECT_ROOT / "data" / "processed"
