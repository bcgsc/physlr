#!/usr/bin/env python3
import igraph as ig
import sys

def read_tsv(g,filename):
	edges = []
	edge_weights = []
	with open(filename) as fin:
		line = fin.readline()
		if line not in ["U\tn\n", "U\tn\tm\n"]:
			print("Unexpected header:", line, file=sys.stderr)
			exit(1)
		reading_vertices = True
		for line in fin:
			if line == "\n":
				line = fin.readline()
				if line == "U\tV\tn\n":
					reading_vertices = False
				else:
					print("Unexpected header:", line, file=sys.stderr)
					exit(1)
				line = fin.readline()
			xs = line.split()
			if reading_vertices:
				if len(xs) == 2:
					g.add_vertex(xs[0], n=int(xs[1]))
				elif len(xs) == 3:
					g.add_vertex(xs[0], n=int(xs[1]), m=int(xs[2]))
				else:
					print("Unexpected row:", line, file=sys.stderr)
					exit(1)
			else:
				if len(xs) == 3:
					edges.append((xs[0], xs[1]))
					edge_weights.append(xs[2])
				else:
					print("Unexpected row:", line, file=sys.stderr)
					exit(1)
	g.add_edges(edges)
	g.es["n"] = edge_weights
	return g

def write_tsv(g,fout):
	if 'm' in g.vs.attributes():
		print("U\tn\tm", file=fout)
	else:
		print("U\tn", file=fout)
	for v in g.vs:
		print(v['name'],end='\t', file=fout)
		for attr in v.attributes():
			if attr != "name":
				print(v[attr],end='\t', file=fout)
		print('', file=fout)
	print("\nU\tV\tn", file=fout)
	for e in g.es:
		l=sorted(list([g.vs[e.source]['name'],g.vs[e.target]['name']]))
		print(*l,sep='\t', end='\t',file=fout)
		for attr in e.attributes():
			print(e[attr], end='\t', file=fout)
		print('', file=fout)

#this function is for when reading graphviz files
"""def sort_vertices(g):
	gsorted=ig.Graph()
	gsorted.add_vertices(sorted(g.vs['name']))
	for e in (g.es):
		gsorted.add_edge(e.source, e.target, n=e['n'])
	return gsorted 
"""

def read_graph(filenames):
	#read_gv = False
	g= ig.Graph()
	for filename in filenames:
		g = read_tsv(g,filename)
	"""
	if read_gv:
		g=sort_vertices(g)
	"""
	return g

def remove_singletons(g):
	singletons=g.vs.select(_degree=0)
	g.delete_vertices(singletons)
	return len(singletons)

def diameter_of_tree(g, weight=None):
	paths=[]
	u=next(iter(g.vs))
	paths=g.get_shortest_paths(u,to=None,weights=weight, mode=ig.ALL,output="vpath")
	path=max(paths,key=len)
	u=path[len(path)-1]
	paths=g.get_shortest_paths(u,to=None,weights=weight, mode=ig.ALL,output="vpath") 
	path=max(paths,key=len)
	v=path[len(path)-1]
	diameter=len(path)-1
	return(u,v,diameter)

def determine_backbones_of_trees(g):
	paths=[]
	for component in g.components():
		gcomponent=g.subgraph(component)
		u,v,_=diameter_of_tree(gcomponent, weight=list(map(int,gcomponent.es['n'])))
		path=gcomponent.get_shortest_paths(u,to=v,weights=list(map(int,gcomponent.es['n'])),mode=ig.ALL,output="vpath")
		paths.append(path)
	paths.sort(key=len, reverse=True)
	return paths

def main():
	filenames=[]
	filenames.extend(sys.argv[1:])
	g = ig.Graph()
	g = read_graph(filenames) #add all filenames into 1 graph
	remove_singletons(g)
	paths= determine_backbones_of_trees(g)
	print(paths)
	fout=open("out",'w')
	write_tsv(g,fout)

if __name__ == "__main__":
	main()
