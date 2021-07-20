import nl4py
import pandas as pd
import multiprocessing
import sys

#print('\n1) Starting the NetLogoControllerServer with: nl4py.startServer()\n')
nl4py.initialize(sys.argv[2])
#print('\n2) Starting the model runs... ')
try:
	n = int(sys.argv[1])
except:
	print('ERROR: Please provide the number of concurrent model runs required as a commandline argument.')

MODEL_PATH = './test14-pantoea.nlogo'
REPORTERS = ['ticks','count turtles']


def init_simulation() -> None:
	global workspace
	workspace = nl4py.create_headless_workspace()
	workspace.open_model(MODEL_PATH)

def run_simulation(run_id : int):
	global workspace
	workspace.command('set efficiency 0.37')
	#workspace.command('set Min/steptime 2')
	workspace.command('set energy_maintenance_pa 0.0015')
	workspace.command('set rep_pa 20')
	workspace.command('set max-time-viability_pa 83')
	workspace.command('set Amonium 18.7')
	workspace.command('set total-length-world 605')
	workspace.command('set colors-gen True')
	workspace.command('set file-id ' + str(run_id))
	workspace.command('setup')
	result = workspace.schedule_reporters(REPORTERS,stop_at_tick=1920)
	return result

if __name__=='__main__':
	print(f'\nOpening {n} copies of the model at {MODEL_PATH} on {multiprocessing.cpu_count()} NetLogo HeadlessWorkspaces')
	names = list(range(n))
	results = []
	with multiprocessing.Pool(processes = multiprocessing.cpu_count(), initializer=init_simulation) as pool:
		results = pd.DataFrame()
		for idx, result in enumerate(pool.map(run_simulation, names)):
			result = pd.DataFrame(result, columns=REPORTERS)
			result['Run'] = idx
			results = results.append(result,ignore_index=True)
	print(f'There are {results.Run.unique().size} runs of results:')
	print(results)