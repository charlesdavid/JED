from pymol import cmd
from pymol.cgo import *
bg_color white
from glob import glob
cd /workspace/cflcyd/git/JED/test/Multi/JED_RESULTS_TEST_MULTI/VIZ/COV
filelist = glob ( "ss_100_Mode_4_COV_frame*.pdb" )
for file in filelist: cmd.load( file, "Mode_4" )
hide lines, Mode_4
show cartoon, Mode_4
cmd.spectrum("b",selection="((all)&*/ca)",quiet=0)
orient
set two_sided_lighting,1
set cartoon_fancy_helices,1
set cartoon_highlight_color,grey50
rebuild
