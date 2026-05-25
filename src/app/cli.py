import os
import yaml
import psutil
import logging
import argparse

from pathlib import Path

from .slurm import submit_sbatch

def skater_annotate(args):
	from annotate import annotate
	if args.config != None:
		with open(args.config,'r') as file:
			config = yaml.safe_load(file)
		annotate(config['annotation']['refFlat'],config['files']['annotation'],config['annotation']['cps_path'],config['files']['RNA-seq'],config['annotation']['use_annotated_introns'],config['annotation']['merge_cps'],config['annotation']['min_reads'],config['annotation']['min_length'],config['annotation']['min_intron_fraction'])
	else:
		annotate(args.refFlat,args.output,args.cps,args.intron,args.use_annotated_introns,args.merge_cps,args.min_intron_reads,args.min_intron_distance,args.min_intron_fraction)

def skater_unpack(args):
	from loadSkater.genome import dumpGenome
	if args.config != None:
		with open(args.config,'r') as file:
			config = yaml.safe_load(file)
		out_path = Path(config['output']['directory'])/'temp/genome.db'
		annotation_path = Path(config['files']['annotation'])
	else:
		annotation_path = Path(args.annotation)
		out_path = Path(args.output)

	dumpGenome(annotation_path,out_path)

def skater_watcher(args):
	# Find path to config file
	config_path = Path(args.config)

	if args.method == 'complete': 
		from skaterWatcher import task_watcher
		task_watcher('skater',config_path)
	elif args.method == 'cleanup':
		from skaterWatcher import task_cleanup
		task_cleanup('skater',config_path)
	else: raise TypeError('Invalid method')

def skater_run(args):
	from optimizeSkater import optimizeGene
	solution = optimizeGene(Path(args.config),args.gene,args.run)
	print(f"Ran\t{solution.scores}")

def skater_load(args):
	from loadSkater import saveGene, LoadingError
	logger = logging.getLogger('__main__')
	try:
		required_events,optional_events = saveGene(args.gene,Path(args.config),False)
		event_str = 'True\t' + ''.join([f'{x},' for x in required_events])[:-1]+';'+''.join([f'{x},' for x in optional_events])[:-1]
		print(event_str)
	except LoadingError as e:
		logger.debug('unsolvable')
		print(e.args[0])

def skater_compile(args):
	from equations.compile import compile

	config = None
	if args.config:
		with open(args.config,'r') as file:
			config = yaml.safe_load(file)
	write_event = config['files']['write_events'] if config else False
	if args.write_event:
		write_event = True

	if args.function_directory:
		function_dir = Path(args.function_directory)
	elif config:
		function_dir = Path(config['files']['function_directory'])
	else: raise TypeError('Need to provide function directory or config path')
	
	print(compile(args.event,function_dir,write_event,False))

def skater_output(args):
	from formats import Bounds
	from output import outputDataframe

	# Define bounds
	bounds = Bounds()
	# gs
	if args.Gs_lb != None: lb = args.Gs_lb
	else: lb = bounds.gs[0]*1.01
	if args.Gs_ub != None: ub = args.Gs_ub
	else: ub = bounds.gs[1]*0.99
	bounds.gs = (lb,ub)
	# gb
	if args.Gb_lb != None: lb = args.Gb_lb
	else: lb = bounds.gb[0]*1.01
	if args.Gb_ub != None: ub = args.Gb_ub
	else: ub = bounds.gb[1]*0.99
	bounds.gb = (lb,ub)
	# spawn
	if args.spawn_lb != None: lb = args.spawn_lb
	else: lb = bounds.spawn[0]*1.01
	if args.spawn_ub != None: ub = args.spawn_ub
	else: ub = bounds.spawn[1]*0.99
	bounds.spawn = (lb,ub)
	# drb
	if args.drb_lb != None: lb = args.drb_lb
	else: lb = bounds.drb[0]*1.01
	if args.drb_ub != None: ub = args.drb_ub
	else: ub = bounds.drb[1]*0.99
	bounds.drb = (lb,ub)
	# splice
	if args.splice_lb != None: lb = args.splice_lb
	else: lb = bounds.splice[0]*1.01
	if args.splice_ub != None: ub = args.splice_ub
	else: ub = bounds.splice[1]*0.99
	bounds.splice = (lb,ub)
	# cleave
	if args.cleave_lb != None: lb = args.cleave_lb
	else: lb = bounds.cleave[0]*1.01
	if args.cleave_ub != None: ub = args.cleave_ub
	else: ub = bounds.cleave[1]*0.99
	bounds.cleave = (lb,ub)
	# elongation
	if args.elongation_lb != None: lb = args.elongation_lb
	else: lb = bounds.elongation[0]*1.01
	if args.elongation_ub != None: ub = args.elongation_ub
	else: ub = bounds.elongation[1]*0.99
	bounds.elongation = (lb,ub)
	# contamination
	if args.contamination_lb != None: lb = args.contamination_lb
	else: lb = bounds.contamination[0]*1.01
	if args.contamination_ub != None: ub = args.contamination_ub
	else: ub = bounds.contamination[1]*0.99
	bounds.contamination = (lb,ub)

	outputDataframe(Path(args.config),Path(args.out),bounds)

class ResourceUsageFilter(logging.Filter):
	def filter(self, record):
		# Add memory usage to the log record
		process = psutil.Process(os.getpid())
		record.mem = int(process.memory_info().rss / (1024 * 1024))  # MB
		record.usr = process.cpu_times().user
		record.sys = process.cpu_times().system
		return True

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
	subparser_watcher.set_defaults(func=skater_watcher)
	subparser_watcher.add_argument('--method','-m',default='complete',help='version of optimizer to use (complete or cleanup)')
	subparser_watcher.add_argument('--config','-c',help='path to config file')

	subparser_load = subparser.add_parser('load',help='loads gene coverage from bam file')
	subparser_load.set_defaults(func=skater_load)
	subparser_load.add_argument('--gene','-g',help = 'gene name')
	subparser_load.add_argument('--config','-c',help = 'path to config file')

	subparser_compile = subparser.add_parser('compile',help='writes and compiles event functions')
	subparser_compile.set_defaults(func=skater_compile)
	subparser_compile.add_argument('--event',help='name of event to be compiled')
	subparser_compile.add_argument('--config','-c',required=False,help = 'path to config file')
	subparser_compile.add_argument('--function_directory',required=False,help='path to function directory') #NOTE: EXPLAIN THIS IS AN OVERIDE OF CONFIG FILE
	subparser_compile.add_argument('--write_event',action='store_true',help='write equations for new events') #NOTE: NEED HELP INFO HERE

	subparser_run = subparser.add_parser('run',help='launch optimization algorithm for a single gene')
	subparser_run.set_defaults(func=skater_run)
	subparser_run.add_argument('--config','-c',help = 'path to config file')
	subparser_run.add_argument('--gene','-g', help = 'gene name')
	subparser_run.add_argument('--run',type=int,help = 'optimization run number')
	
	subparser_output = subparser.add_parser('output',help='save optimizer output')
	subparser_output.set_defaults(func=skater_output)
	subparser_output.add_argument('--out','-o',help='path to save output file')
	subparser_output.add_argument('--config','-c',help='path to config file')
	subparser_output.add_argument('--Gs_lb',required=False,type=float,help='lower bound for Gs')
	subparser_output.add_argument('--Gs_ub',required=False,type=float,help='upper bound for Gs')
	subparser_output.add_argument('--Gb_lb',required=False,type=float,help='lower bound for Gb')
	subparser_output.add_argument('--Gb_ub',required=False,type=float,help='upper bound for Gb')
	subparser_output.add_argument('--spawn_lb',required=False,type=float,help='lower bound for spawn rate')
	subparser_output.add_argument('--spawn_ub',required=False,type=float,help='upper bound for spawn rate')	
	subparser_output.add_argument('--drb_lb',required=False,type=float,help='lower bound for drb release rate')
	subparser_output.add_argument('--drb_ub',required=False,type=float,help='upper bound for drb release rate')
	subparser_output.add_argument('--splice_lb',required=False,type=float,help='lower bound for splicing rates')
	subparser_output.add_argument('--splice_ub',required=False,type=float,help='upper bound for splicing rates')
	subparser_output.add_argument('--cleave_lb',required=False,type=float,help='lower bound for cleavage rate')
	subparser_output.add_argument('--cleave_ub',required=False,type=float,help='upper bound for cleavage rate')
	subparser_output.add_argument('--elongation_lb',required=False,type=float,help='lower bound for elongation rate')
	subparser_output.add_argument('--elongation_ub',required=False,type=float,help='upper bound for elongation rate')
	subparser_output.add_argument('--contamination_lb',required=False,type=float,help='lower bound for contamination')
	subparser_output.add_argument('--contamination_ub',required=False,type=float,help='upper bound for contamination')

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
		cmd = f"skater {temp['command']}"
		
		for flag,val in temp.items():
			if flag in ['slurm','p','partition','job_name','time_limit','time','mem','cpu','log','func','command']: 
				# Skip commands related to slurm submission
				continue
			elif isinstance(val,bool): 
				if val: cmd += f" --{flag}"
			elif val == '' or val == None: continue
			elif len(flag) == 1: cmd += f" -{flag} {val}"
			else: cmd += f" --{flag} {val}"
		submit_sbatch(args.partition,args.job_name,args.mem,args.time_limit,int(args.cpu),Path(args.log),cmd)
		print(cmd)
	else: args.func(args)
	
if __name__ == '__main__':
	main()