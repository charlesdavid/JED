from pymol import cmd
from pymol.cgo import *
bg_color white
from glob import glob
cd /workspace/cflcyd/git/JED/test/Single/JED_RESULTS_TEST_SINGLE/VIZ_TEST6
filelist = glob ( "ss_151_Essential_Modes_1_10_CORR_frame*.pdb" )
for file in filelist: cmd.load( file, "Essential_Modes_10" )
hide lines, Essential_Modes_10
show cartoon, Essential_Modes_10
cmd.spectrum("b",selection="((all)&*/ca)",quiet=0)
orient
set two_sided_lighting,1
set cartoon_fancy_helices,1
set cartoon_highlight_color,grey50
rebuild
