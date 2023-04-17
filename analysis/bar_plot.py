from matplotlib.pyplot import savefig, subplots, show, title, suptitle
from numpy import array, zeros, arange

nodes_number = 4

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

lines = 0  
section_witdh = nodes_number * 2
data = []

with open('data/times.dat') as file:	
	for line in file.readlines()[1:]:
		data += [[num(x) for x in line.split()[1:]]]
		lines += 1

#data = data[1:]
sections = lines // section_witdh

for s in range(sections):
	mpi = data[s * section_witdh : nodes_number + s * section_witdh]
	acc = data[nodes_number + s * section_witdh : 2 * nodes_number + s * section_witdh]
	dim = mpi[0][0]
	itr = mpi[0][1]
	versions = (mpi, acc)
	computation = [array([t for *_, t, _ in version]) for version in versions]
	communication = [array([t for *_, t in version]) for version in versions]

	measures = (communication, computation)
	nodes = array([2 ** i for i in range(nodes_number)])

	fig, ax = subplots(layout='constrained')
	bottom = zeros(nodes_number)
	width = 0.2
	multiplier = 0
	i = 0
	bar_colors = ('black', 'black', 'red', 'yellow')
	label = ("Communication", None, "MPI computation", "OpenACC+MPI computation")

	#print(measures)
	for measure_type in measures:
	  
		for times in measure_type:
			if i > 1:
				bottom = measures[0][i - 2]
		  
			offset = width * multiplier - width

			p = ax.bar(nodes + offset, times, width, label = label[i], bottom = bottom, \
						     color = bar_colors[i], edgecolor = 'black')
						     
			i += 1
			multiplier += 1
		
		multiplier = 0
		bottom = zeros(nodes_number)

	ax.set_xticks(nodes)
	ax.set_ylabel(r"Time [$s$]")
	ax.set_xlabel("Number of nodes")
	ax.legend(loc = "upper right")
	ax.grid(linestyle = '--', axis = 'y')

	title(f"Size of grid: {dim}"+r"$\times$"+f"{dim}\nIterations: {itr}", fontsize = 10)
	suptitle('Communication and computation times per number of nodes', fontsize = 13, y = 1.03, x = 0.54)
	savefig(f'analysis/{dim}_{itr}.png')
    # show()

