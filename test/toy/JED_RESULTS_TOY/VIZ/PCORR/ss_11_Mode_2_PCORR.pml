from pymol import cmd
from pymol.cgo import *
bg_color white
from glob import glob
filelist = glob ( "ss_11_Mode_2_PCORR_frame*.pdb" )
for file in filelist: cmd.load( file, "Mode_2" )
hide lines, Mode_2
show cartoon, Mode_2
cmd.spectrum("b",selection="((all)&*/ca)",quiet=0)
orient
set two_sided_lighting,1
set cartoon_fancy_helices,1
set cartoon_highlight_color,grey50
util.performance(0)
rebuild