import nl4py
import pandas as pd
import sys

print("\n1) Starting the NetLogoControllerServer with: nl4py.startServer()\n")
nl4py.startServer(sys.argv[2])
print('\n2) Starting the model runs... ')
try:
	n = sys.argv[1]
except:
	print("ERROR: Please provide the number of concurrent model runs required as a commandline argument.")

model = "./test14-pantoea.nlogo"
print('\n2.1) Creating ' + n + ' NetLogo HeadlessWorkspaces with: workspace = nl4py.newNetLogoHeadlessWorkspace()')
print('\n2.2) Opening ' + str(n) + ' copies of the model at ' + model + ' on the NetLogo HeadlessWorkspaces with: workspace.openModel("model")')
for i in range(0,int(n)):
	workspace = nl4py.newNetLogoHeadlessWorkspace()
	workspace.openModel(model)

print("\n2.3) Get all workspaces back with: workspaces = nl4py.getAllExistingWorkspaces() \n\tSetting the parameters for all " + str(n) + " models to random values with workspace.setParamsRandom()")
for idx, workspace in enumerate(nl4py.getAllHeadlessWorkspaces()):
	workspace.setParamsRandom()
	workspace.command("set efficiency 0.37")
#	workspace.command("set Min/steptime 2")
	workspace.command("set energy_maintenance_pa 0.0015")
	workspace.command("set rep_pa 20")
	workspace.command("set max-time-viability_pa 83")
	workspace.command("set Amonium 18.7")
	workspace.command("set total-length-world 605")
	workspace.command("set colors-gen True")
	workspace.command("set file-id " + str(idx))

print('\n2.4) Send setup and go commands to each model using: workspace.command("setup") and workspace.command("go") ')
results = pd.DataFrame()
reporters = ["ticks"]
#reporters = ["ticks","count turtles"]
for ir in range(0,6):
	for idx, workspace in enumerate(nl4py.getAllHeadlessWorkspaces()):
		workspace.command("set file-ir " + str(ir))
		workspace.command("setup")
		workspace.scheduleReportersAndRun(reporters,stopAtTick=1920)
		print("dosp, idx, ir ", idx,ir)
	for idx, workspace in enumerate(nl4py.getAllHeadlessWorkspaces()):
		result = workspace.awaitScheduledReporterResults()
		result = pd.DataFrame(result, columns=reporters)
		result["workspace"] = idx
		result["iteration"] = ir
		results = results.append(result)

results = results.set_index(["workspace","iteration"])
print(results)
print('\n3) Shutdown the server to release compute resources using: nl4py.stopServer()')
nl4py.stopServer()
print('\n\n------------------------ Thanks for trying NL4PY -------------------------\n')
