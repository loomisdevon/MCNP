# Dual Rotating Radiation Scattering Mask MCNP Code
# Rotate neutron and gamma flux from 0 to 2PI around mask and return signal received from detector 

# Authors: Ivan Novikov and Devon Loomis



import os
import os.path
import math
import fileinput
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import time
import math
from decimal import Decimal
import csv
from tqdm import *


###################################### creating all MCNP input files for simulation of the increasing relative distance between source and detector #######################
''' 														Parameters:
														file_name: MCNP input template
														y-pos: y position of center of circle of source path
														r: radius of circle of source path
														init: inital angle away from xy-plane
														final angle away from xy-plane
														angle step size
														limit: longest distance between the two
'''

def createFiles(file_name, z_pos, r, init, final, step_size):
	fileList = []
	marker=0
	for new_ang in range(init, final, step_size):

		text_search = None
		f =open(file_name)
		for line in f:
			words = line
			sdef = words[0:4]
			if (sdef == "SDEF"):
				text_search = words
				break
		f.close()
		
		rad_ang = math.radians(new_ang)

		x_pos = round(r * np.cos(rad_ang),5)
		y_pos = round(r * np.sin(rad_ang),5)
		vecx_pos = round(-1 * np.cos(rad_ang),5)
		vecy_pos = round(-1 * np.sin(rad_ang),5)
		theta_rad = np.arctan(z_pos/r)
		vecz_pos = round(-1 * (theta_rad/(np.pi/2)),5)
		replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 par=n" + "\n"
		read_name = file_name
		write_name = "inp" + str(new_ang) + ".txt"

		f1 = open(read_name, 'r')
		f2 = open(write_name, 'w')

		for lines in f1:
			f2.write(lines.replace(text_search, replacement_text))

		f1.close()
		f2.close()
		fileList.append(write_name)

	return (fileList)
		


################################# delete runtpe files after every set of commands and delete all output files and input files after program run #######################
'''																Parameters
															directory: directory containing all files 
															file: KSEF_2 #####################
															remove_all: test to determine whether to delete all files or only runtpe files
'''

def removeFiles(directory, file, outfile, remove_all):
	dir_name = directory
	for fname in os.listdir(dir_name):
		if fname.startswith("binRun"):
			os.remove(os.path.join(dir_name, fname))
		if (fname.startswith(file[:-4]) or fname.startswith(outfile[:-4])) and remove_all:
			if (fname != file):
				os.remove(os.path.join(dir_name, fname))



####################### read MCNP output file, find and return flux value #########################
#######################_file_: MCNP output file name ##################################

def readFlux(_file_):
	flux_ = 0
	error_ = 0
	with open(_file_, 'r') as outfile:
		for line in outfile:
			if ('+                                   *Neutron flux in detector*' in line):
				lines = [outfile.readline() for i in range(9)] #this reads 9 lines after the fc4 comment
				spectrum = [outfile.readline() for i in range(13)] #this reads 13 lines which contain spectrum
	            #each line has an index [0]-[12]
	    		#print(type(spectrum[1]))
	    		#print(spectrum[1])

	    		#print(float(spectrum[1].split()[1])) #this splits spectrum[i] using spaces
	            #each spectrum[i].split() has three new indeces [0]-[2]
	            #float converts each string to float
	            #Neutron energy is in [0]
	            #Neutron counts are in [1]
	            #Error is in [2]
				#tmp = 0.0
				#print (spectrum)
				for j in range(13):
					flux_ += float(spectrum[j].split()[1])
					error_ += float(spectrum[j].split()[2])

            #Fluxin3[i] = tmp

	return flux_, error_


#**********************MAIN**************************
dir_ = 'C:\\Users\\devon\\Documents\\DRRS Mask\\MCNP\\'
file_ = 'inp.txt'
outFile_ = "out.txt"
file_name_ = dir_ + file_
outFile_name_ = dir_ + outFile_




packet = 16 # no. of simultaneous commands sent to MCNP shell


#start = time.time()

activity = 3.7 * (10**10)
nps = 1 * (10**7)
radius = 45 # radius of circle of source path
init_ang = 0 # inital angle away from xy-plane
final_ang = 360 # final angle away from xy-plane
step = 2 # angle step siz
intensity = 1*10**7
t = 100
z = [0,10]

originalZCountsArray = []
originalZErrorArray = []


for z_ in z:
	start = time.time()
	removeFiles(dir_, file_, outFile_, True) # purge directory of any existing MCNP files from previous run
	files = createFiles(file_name_, z_, radius, init_ang, final_ang, step) # create all MCNP input files

	commands = []
	outFileList = []
	j = init_ang

	#create set of commands for subprocess of all input files
	for i in range(int((final_ang - init_ang) / step)):
		binFile = "binRun" + str(j) + ".r"	
		outFile = ("out" + str(j) + ".txt")
		commands.append("mcnp6 i=" + files[i] + " o=" +  outFile + " runtpe=" + binFile)
		outFileList.append(outFile)
		j += step



	print("Simulating...")

	# give subprocess pak amount of parallel programs to execute until all commands are executed 
	for x in tqdm(range(0,int((final_ang - init_ang) / step),(packet))):	
		if (x < (len(commands) - packet)):
			commandsub = commands[x:(x+packet)]
		else:
			commandsub = commands[x:]
		processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, cwd=dir_) for cmd in commandsub]
		removeFiles(dir_, file_, outFile_, False) # remove runtpe files

		for p in processes:
			p.wait()

	ang = init_ang
	
	fluxList = []
	errorList = []
	sourceAngleList = []


	############################read and gather flux values and source distances for each output file and add them to lists###################################
	for f in outFileList:
		flux, error = readFlux(f)

		fluxList.append(flux)
		errorList.append(error)
		rad_ang = math.radians(ang)
		sourceAngleList.append(rad_ang)
		ang += step

	#removeFiles(dir_, file_, outFile_,True)
	end = time.time()
	print("Runtime: ", round((end - start)/60, 2), " mins")
	
	

	fluxArray = np.array(fluxList)
	errorArray = np.array(errorList)
	angleArray = np.array(sourceAngleList)
	

	countsArray = fluxArray * intensity * t


	countsSum = np.sum(countsArray)


	normalizedCountsArray = countsArray / countsSum

	deltaNorm = np.sqrt(np.sum((countsArray*errorArray)**2))

	deltaN_N2 = (normalizedCountsArray**2)*(((errorArray)**2) + ((deltaNorm / countsSum)**2))

	deltaN_Stat_N2 = (normalizedCountsArray**2)*((1/np.sqrt(countsArray))+((np.sqrt(countsSum)/countsSum)**2))

	normalizedCountsErrorArray = np.sqrt(deltaN_N2 + deltaN_Stat_N2)


	


	if z_ == z[0]:
		originalZCountsArray = normalizedCountsArray
		originalZErrorArray = normalizedCountsErrorArray
		DeltaN = 0
		DeltaNError = 0
	else:
		DeltaNArray = (normalizedCountsArray - originalZCountsArray)**2
		DeltaN = np.sum(DeltaNArray)
		DeltaNError = np.sqrt(np.sum(((2*DeltaNArray*np.sqrt(normalizedCountsErrorArray**2 + originalZErrorArray**2))/(np.abs(normalizedCountsArray - originalZCountsArray)))**2))






	titlestring = "z = " + str(z_)
	title = [titlestring]
	headers = ["Angle","Raw Flux","Raw Error", "Normalized Counts", "Normalized Counts Error", "MCNPERROR", "STATISTICALERROR"]
	deltaResults = [DeltaN, DeltaNError]
	totalArray= [angleArray, fluxArray, errorArray, normalizedCountsArray, normalizedCountsErrorArray, deltaN_N2, deltaN_Stat_N2]

	with open("data_output_config1VarianceRed_010.csv", "a", newline='') as out:
		writer = csv.writer(out)
		writer.writerow(title)
		writer.writerow(headers)
		writer.writerows(zip(*totalArray))
		writer.writerow(deltaResults)
	


	print("DeltaN = " + str(round(DeltaN,8)) + " DeltaNError = " + str(round(DeltaNError,8)))
	plt.plot(sourceAngleList, normalizedCountsArray, label ='S(z=' + str(z_) + ")")
	#plt.errorbar(sourceAngleList, normalizedCountsArray, normalizedCountsErrorArray, label ='S(z=' + str(z_) + ")")

################################Plot###################################################
plt.xlabel('phi (rad)')
plt.ylabel('neutron/cm^2/source neutron')
#plt.title('All-Poly (M1) configuration     \u0394'+'N = ' + str(f"{Decimal(DeltaN):.3E}") + '     \u03C3'+' = ' + str(f"{Decimal(sigmaN):.3E}"))
plt.autoscale(enable=True, axis='both',tight=False)
plt.legend(loc='upper right')

plt.savefig('configuration1_countsVarianceRed_010_withoutput.png')
plt.show()



###########################END MAIN###############################