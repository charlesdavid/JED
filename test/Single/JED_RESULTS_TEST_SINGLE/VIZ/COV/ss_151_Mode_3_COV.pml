from pymol import cmd
from pymol.cgo import *
bg_color white
from glob import glob
cd /workspace/cflcyd/git/JED/test/Single/JED_RESULTS_TEST_SINGLE/VIZ/COV
filelist = glob ( "ss_151_Mode_3_COV_frame*.pdb" )
for file in filelist: cmd.load( file, "Mode_3" )
hide lines, Mode_3
show cartoon, Mode_3
cmd.spectrum("b",selection="((all)&*/ca)",quiet=0)
orient
set two_sided_lighting,1
set cartoon_fancy_helices,1
set cartoon_highlight_color,grey50
rebuild
