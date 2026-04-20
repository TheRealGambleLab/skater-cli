import yaml
import logging
import argparse
import numpy as np

from pathlib import Path

from .slurm import submit_sbatch
from .cli import skater_annotate, skater_unpack, ResourceUsageFilter

def tt_watcher(args):
	# Find path to config file
	config_path = Path(args.config)

	if args.method == 'complete': 
		from skaterWatcher import task_watcher
		task_watcher('tt-skater',config_path)
	elif args.method == 'cleanup':
		from skaterWatcher import task_cleanup
		task_cleanup('tt-skater',config_path)
	else: raise TypeError('Invalid method')

def tt_run(args):
	from optimizeTT import optimizeGene
	solution = optimizeGene(Path(args.config),args.gene,args.run)
	print(f"Ran\t{np.sum(solution.scores)}")

def tt_load(args):
	from loadSkater import saveGene, LoadingError
	logger = logging.getLogger('__name__')
	try:
		required_events,optional_events = saveGene(args.gene,Path(args.config),True)
		event_str = 'True\t' + ''.join([f'{x},' for x in required_events])[:-1]+';'+''.join([f'{x},' for x in optional_events])[:-1]
		print(event_str)
	except LoadingError as e:
		logger.debug('unsolvable')
		print(e.args[0])

def tt_compile(args):
	from equations.compile import compile

	write_event = args.write_event
	if args.function_directory:
		function_dir = Path(args.function_directory)
	elif args.config:
		with open(args.config,'r') as file:
			config = yaml.safe_load(file)
		function_dir = Path(config['files']['function_directory'])
		write_event = config['files']['write_events']
	else: raise TypeError('Need to provide function directory or config path')
	
	print(compile(args.event,function_dir,write_event,True))

def tt_output(args):
	from output import outputDataframe
	outputDataframe(Path(args.config),Path(args.out))

def tt_pro(args):
	from loadSkater.pro import getProCov
	# Get config file
	with open(args.config,'r') as file:
		config:dict = yaml.safe_load(file)

	getProCov(config,Path(args.out))

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--resource_monitor',action='store_true')
	parser.add_argument('--slurm',action='store_true')
	parser.add_argument('--partition','-p',default='unlimited')
	parser.add_argument('--job_name',default='skater')
	parser.add_argument('--time_limit','--time',default='72:00:00',type=str)
	parser.add_argument('--mem',default='20gb')
	parser.add_argument('--cpu',default=4,type=int)
	parser.add_argument('--log',default='./slurm-%j.out')

	subparser = parser.add_subparsers(dest='command')
	
	subparser_annotate = subparser.add_parser('annotate',help='creates gene annotation file')
	subparser_annotate.set_defaults(func=skater_annotate)
	subparser_annotate.add_argument('--refFlat',default=None,help='path to annotation refflat')
	subparser_annotate.add_argument('--config',default=None,help='path to config file')
	subparser_annotate.add_argument('--output',default=None,help='path to output tsv')
	subparser_annotate.add_argument('--cps',default=None,help='path to CPS annotation file')
	subparser_annotate.add_argument('--intron',default=None,help='path to bed file for intron definition')
	subparser_annotate.add_argument('--use_annotated_introns',action='store_true',help='use refflat intron annotations')
	subparser_annotate.add_argument('--merge_cps',action='store_true',help='keep both refflat and cps file annotations')
	subparser_annotate.add_argument('--min_intron_reads',default=2,type=int,help='minimum number of reads required to support an intron definition')
	subparser_annotate.add_argument('--min_intron_distance',default=50,type=int,help='minimum distance required to define an intron')
	subparser_annotate.add_argument('--min_intron_fraction',default=0.05,type=float,help='minimum fraction of total junctions to define intron')
	# TODO: ADD HELP LINE TO SUMMARIZE COLUMN NAMES OF OUTPUT

	subparser_unpack = subparser.add_parser('unpack',help='unpacks data in annotation file and stores in pickled database')
	subparser_unpack.set_defaults(func=skater_unpack)
	subparser_unpack.add_argument('--annotation','-a',default=None,help='path to annotation tsv')
	subparser_unpack.add_argument('--output','-o',default=None,help='path to output database file')
	subparser_unpack.add_argument('--config','-c',default=None,help='path to config file')

	subparser_watcher = subparser.add_parser('watcher',help = 'launch watcher script to submit jobs to run on HPC')
	subparser_watcher.set_defaults(func=tt_watcher)
	subparser_watcher.add_argument('--method','-m',default='complete',help='version of optimizer to use (complete or cleanup)')
	subparser_watcher.add_argument('--config','-c',help='path to config file')

	subparser_load = subparser.add_parser('load',help='loads gene coverage from bam file')
	subparser_load.set_defaults(func=tt_load)
	subparser_load.add_argument('--gene','-g',help = 'gene name')
	subparser_load.add_argument('--config','-c',help = 'path to config file')

	subparser_compile = subparser.add_parser('compile',help='writes and compiles event functions')
	subparser_compile.set_defaults(func=tt_compile)
	subparser_compile.add_argument('--event',help='name of event to be compiled')
	subparser_compile.add_argument('--config','-c',required=False,help = 'path to config file')
	subparser_compile.add_argument('--function_directory',required=False,help='path to function directory') #NOTE: EXPLAIN THIS IS AN OVERIDE OF CONFIG FILE
	subparser_compile.add_argument('--write_event',action='store_true',help='write equations for new events') #NOTE: NEED HELP INFO HERE

	subparser_run = subparser.add_parser('run',help='launch optimization algorithm for a single gene')
	subparser_run.set_defaults(func=tt_run)
	subparser_run.add_argument('--config','-c',help = 'path to config file')
	subparser_run.add_argument('--gene','-g', help = 'gene name')
	subparser_run.add_argument('--run',type=int,help = 'optimization run number')
	
	subparser_output = subparser.add_parser('output',help='save optimizer output')
	subparser_output.set_defaults(func=tt_output)
	subparser_output.add_argument('--out','-o',help='path to save output file')
	subparser_output.add_argument('--config','-c',help='path to config file')

	subparser_pro = subparser.add_parser('pro',help='save PROseq coverage to file')
	subparser_pro.set_defaults(func=tt_pro)
	subparser_pro.add_argument('--out','-o',help='path to save output file (.bw are split into pos and negative strand)')
	subparser_pro.add_argument('--config','-c',help='path to config file')
	args = parser.parse_args()

	# Set up logger
	logger = logging.getLogger('__main__')
	hdr = logging.StreamHandler()
	if args.resource_monitor:
		fmt = logging.Formatter('%(asctime)s\t%(message)s\t %(mem)s MB\t %(usr)s sec %(sys)s sec','%Y-%m-%d %H:%M:%S')
		logger.addFilter(ResourceUsageFilter())
	else:
		fmt = logging.Formatter('%(asctime)s\t%(message)s','%Y-%m-%d %H:%M:%S')
	hdr.setFormatter(fmt)
	logger.addHandler(hdr)
	logger.setLevel('INFO')


	if args.command is None: # No command given
		# Print help info
		parser.print_help()
	elif args.slurm: # Submit job as a slurm submission
		# Unpack arguments
		temp = vars(args)

		# Write command
		cmd = f"tt-skater {temp['command']}"
		
		for flag,val in temp.items():
			if flag in ['slurm','p','partition','job_name','time_limit','time','mem','cpu','log','func','command']: 
				# Skip commands related to slurm submission
				continue
			elif isinstance(val,bool): 
				if val: cmd += f" --{flag}"
			elif val == '': continue
			elif len(flag) == 1: cmd += f" -{flag} {val}"
			else: cmd += f" --{flag} {val}"
		submit_sbatch(args.partition,args.job_name,args.mem,args.time_limit,int(args.cpu),Path(args.log),cmd)
		print(cmd)
	else: args.func(args)
	
if __name__ == '__main__':
	main()