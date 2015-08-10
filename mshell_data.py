#!/usr/bin/env
'''
		University of Southampton	August 2015

			mshell_data.py 

mshell_data reads in atomic data from m-shell - the atomic data set from
Stuart Sim's code, and converts this to python format

Usage:
	python mshell_data prefix folder

	prefix gives determines the output filenames, e.g. xxxx_lines.py, xxxx_phot.pt
	folder tells you where the files ATOM.MODELS etc. are stored.

History:
	1508 JM -- Coded
'''
import numpy as np 
import atomic_classes as cls 
import os, sys 
import atomic_sub as sub 
'''
Format of stuarts files

ATOM.MODELS (Level info)
PHIXS (Photoionization cross-sections)
LINELIST (List of lines)
'''


# A few constants we use in our calculations
MAX = 1e8
A21_CONSTANT=7.429297e-22
ANGSTROM = 1e-8
HEV	=4.13620e-15	# Planck's constant in eV 
C =2.997925e10	
MEGABARN = 1e-18
default_folder = "m-shell"

def read_atom_model(fname="%s/ATOM.MODELS" % default_folder, outname="test.out", z_select= None, write=True):

	'''
	read_atom_model reads Stuarts data from the ATOM.MODELS file and converts it into
	an array of cls.level class instances to be written out to file

	Arguments:
		fname		string
					fname to read in, e.g. m-shell/ATOM.MODELS
		outname		string

	''' 

	# read in the data and store in the data array
	f = open(fname, "r")

	nions = 0
	data = []
	for line in f:
		d = line.split()
		if len(d) == 4:
			nions += 1

		if len(d) > 0:
			data.append(d)
	f.close()

	print "READ %i IONS FROM %s" % (nions, fname)
	print "READ %i TOTAL LEVELs FROM %s" % (len(data)-nions, fname)

	i = 0

	# create a blank array to store classes
	levels = []
	
	# cycle through the data and place in class instances
	while i < len(data) and i < MAX:

		if len(data[i]) == 4:		# we have a summary record
			z = int(data[i][0])
			ion = int(data[i][1])
			nlevs = int(data[i][2])
			threshold = float(data[i][3])

			for j in range(nlevs):
				i+=1
				lvl = float(data[i][0])
				energy_ev = float(data[i][1])
				g = float(data[i][2])
				rad_rate = 1e-20
				nstring="nn"
				ionpot = energy_ev - threshold
			
				if z_select == z or z_select==None:
					lev = cls.level (z, ion, lvl, ionpot, energy_ev, g, rad_rate, nstring, 0)
					levels.append(lev)

		i+=1

	levels = np.array(levels)

	# write out the file
	if write:
		sub.write_level_file(levels, outname, levmax = 1e50, append = False, z = None, ion = None)

	return levels


def read_line_list(level_class_array, fname="%s/LINELIST"  % default_folder, outname="lines.out", z_select = None, write=True):

	'''
	read_atom_model reads Stuart's data from the ATOM.MODELS file and converts it into
	an array of cls.level class instances to be written out to file
	''' 

	f = open(fname, "r")

	nions = 0
	data = []
	for line in f:
		d = line.split()
		if len(d) == 4:
			nions += 1

		if len(d) > 0:
			data.append(d)

	f.close()

	print "WE HAVE READ %i lines" % len(data)

	i = 0

	# create a blank array to store classes
	lines = []
	# (l.z, l.ion, l.lvl, l.ionpot, l.E, l.g, l.rad_rate, l.nnstring))
	while i < len(data):

		if len(data[i]) == 3:		# we have a summary record
			z = int(data[i][0])
			ion = int(data[i][1])
			nlines = int(data[i][2])

			for j in range(nlines):
				i+=1
				nline = int(data[i][0])
				lower = int(data[i][1])
				upper = int(data[i][2])
				A_value = float(data[i][3])

				# only bother if we want this element
				if z_select == z or z_select==None:
					found_lower = False
					found_upper = False

					for ilev in range(len(level_class_array)):
						if level_class_array[ilev].z == z and level_class_array[ilev].ion == ion:
							if level_class_array[ilev].lvl == lower:
								El = level_class_array[ilev].E
								gl = level_class_array[ilev].g
								found_lower = True

							if level_class_array[ilev].lvl == upper:
								Eu = level_class_array[ilev].E
								gu = level_class_array[ilev].g
								found_upper = True


					if found_upper * found_lower == False:
						print "NOT FOUND! Line %i %i %i %i" % (z, ion, j, nlines)

					energy_gap = Eu - El
					freq = energy_gap / HEV				# get frequency of line

					try:
						wavelength = C / freq / ANGSTROM	# get wavelength in ANGSTROMS
					except ZeroDivisionError:
						print "ZeroDivisionError: Line %i %i %i %i %8.4e %8.4e" % (z, ion, j, nlines, energy_gap, freq)
						wavelength = 1e5
					#except ZeroDivisionError:
					#wavelength = 0	# get wavelength in ANGSTROMS
					

					# get f value from A value
					try:
						f = A_value / A21_CONSTANT / gl * gu / freq / freq 
					except ZeroDivisionError:
						f = 0.0000

					lines.append(cls.line (z, ion, wavelength, freq, f, gl, gu, El, Eu, lower, upper))

		i+=1

	lines = np.array(lines)

	# write out the file
	if write:
		sub.write_line_file(lines, outname, levmax = 1e50, append = False, z = None, ion = None)


def read_xsections(level_class_array, fname="%s/PHIXS" % default_folder,outname="xs.out", z_select = None, write=True):

	'''
	read_atom_model reads Stuart's photoionization data from the PHIXS file and converts it into
	an array of cls.top_photo_mac class instances to be written out to file.

	Arguments:
		level_class_array	array-like
							array of level class instances linked to data to 
							get thresholds


	''' 

	f = open(fname, "r")

	nxs = 0
	data = []
	for line in f:
		d = line.split()
		if len(d) == 6:
			nxs += 1

		if len(d) > 0:
			data.append(d)

	f.close()

	print "WE HAVE READ %i XSECTIONS" % nxs

	i = 0

	# create a blank array to store classes
	xs = []
	# (l.z, l.ion, l.lvl, l.ionpot, l.E, l.g, l.rad_rate, l.nnstring))
	while i < len(data) and i < MAX:

		if len(data[i]) == 6:		# we have a summary record
			z = int(data[i][0])
			lower_ion = int(data[i][1])
			lower_level = int(data[i][2])
			upper_ion = int(data[i][3])
			upper_level = int(data[i][4])
			entries = int(data[i][5])

			if z == z_select or z_select == None:
				# find the threshold for this ion and level
				for ilev in range(len(level_class_array)):
					if level_class_array[ilev].z == z and level_class_array[ilev].ion == lower_ion:
						if level_class_array[ilev].lvl == lower_level:
							threshold = -level_class_array[ilev].ionpot



				XS = np.zeros(entries)
				energy = np.zeros(entries)
				for j in range(entries):
					i+=1
					energy[j] = threshold + float(data[i][0])
					XS[j] = float(data[i][1]) * MEGABARN

				#print i

				xs.append(cls.top_photo_mac(z, lower_ion, lower_level, upper_level, threshold, entries, energy, XS))

		i+=1

	# write out the file
	if write:
		sub.write_top_macro(xs, outname, suffix = "Mac", append = False, levmax = 100)

# Next lines permit one to run the routine from the command line with various options -- see docstring
if __name__ == "__main__":

	if len(sys.argv) == 3:
		prefix = sys.argv[1]
		folder = sys.argv[2]
	else: 
		print __doc__
		sys.exit()

	levels = read_atom_model(fname="%s/ATOM.MODELS" % folder, outname="%s_levels.py" % prefix)
	read_line_list(levels, fname="%s/LINELIST" % folder, outname="%s_lines.py" % prefix)
	read_xsections(levels, fname="%s/PHIXS" % folder, outname="%s_phot.py" % prefix)


