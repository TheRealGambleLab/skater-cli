import subprocess

from pathlib import Path

def submit_sbatch(partition:str,name:str,mem:str,time:str,cpu:int,out_path:Path,cmd:str) -> subprocess.CompletedProcess:
	batch_command = f"#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH --job-name={name}\n#SBATCH --mem={mem}\n#SBATCH --time={time}\n#SBATCH --cpus-per-task={cpu}\n#SBATCH --open-mode=append\n#SBATCH --output={out_path.resolve()}\n{cmd}"
	proc = subprocess.run(['echo', '-e',batch_command], check=True,text=True,capture_output=True)
	output = subprocess.run(['sbatch'],input=proc.stdout,capture_output=True,text=True)
	return output