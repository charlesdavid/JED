from pymol import cmd
from pymol.cgo import *
bg_color white
from glob import glob
filelist = glob ( "ss_11_Essential_Modes_3_COV_frame*.pdb" )
for file in filelist: cmd.load( file, "Essential_Modes_3" )
hide lines, Essential_Modes_3
show cartoon, Essential_Modes_3
cmd.spectrum("b",selection="((all)&*/ca)",quiet=0)
orient
set two_sided_lighting,1
set cartoon_fancy_helices,1
set cartoon_highlight_color,grey50
util.performance(0)
rebuild
