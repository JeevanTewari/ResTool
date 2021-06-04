from Bio import Phylo

workspace = '/Users/jeevantewari/Desktop/NewProj/OR6Y1_cont/tree'

tree = Phylo.read(workspace, "newick")
tree.ladderize()  # Flip branches so deeper clades are displayed at top
Phylo.draw(tree)