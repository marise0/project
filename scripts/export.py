import sys
import os
sys.path.append(os.path.abspath(os.path.join('..')))
from config import *
print(f"export VCF_PATH={vcf_path}")
print(f"export OUTPUT_PATH={output_path}")
print(f"export SCRIPT_PATH={script_path}")   
print(f"export DOC={doc}")
