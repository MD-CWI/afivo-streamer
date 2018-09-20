# Gnuplot script to make plots of the convergence tests
# Usage: gnuplot plot.gp
# Output: pdf files

# This can be used to plot the stored results, or newly generated data
mypath="output"
# mypath="stored_output"

set style data lines
# set key bottom center
set key outside
set terminal pdf enhanced size 15cm,9cm

set xlabel "t (ns)"
set xrange [1.5:]

set output "cyl_pos_grid.pdf"
set title "Cylindrical (no photoi): effect of grid"
plot mypath."/cyl_pos_adx2.0_log.txt" u ($2*1e9):7 title "adx=2.0",\
     mypath."/cyl_pos_adx1.0_log.txt" u ($2*1e9):7 title "adx=1.0",\
     mypath."/cyl_pos_adx0.5_log.txt" u ($2*1e9):7 title "adx=0.5",\
     mypath."/cyl_pos_adx0.25_log.txt" u ($2*1e9):7 title "adx=0.25",\
     mypath."/cyl_pos_unif_6.0_log.txt" u ($2*1e9):7 dt 2 title "unif. 6.0mu",\
     mypath."/cyl_pos_unif_3.0_log.txt" u ($2*1e9):7 dt 2 title "unif. 3.0mu",\
     mypath."/cyl_pos_unif_1.5_log.txt" u ($2*1e9):7 dt 2 title "unif. 1.5mu"

set output "cyl_vs_3d_pos.pdf"
set title "Cylindrical vs 3D"
plot mypath."/cyl_pos_adx1.0_log.txt" u ($2*1e9):7 title "cyl (no photoi), adx=1.0",\
     mypath."/cyl_pos_adx0.25_log.txt" u ($2*1e9):7 title "cyl (no photoi), adx=0.25",\
     mypath."/3d_pos_adx1.0_log.txt" u ($2*1e9):7 dt 2 title "3d (no photoi), adx=1.0",\
     mypath."/cyl_pos_ph_adx1.0_log.txt" u ($2*1e9):7 title "cyl (photoi), adx=1.0",\
     mypath."/cyl_pos_ph_adx0.25_log.txt" u ($2*1e9):7 title "cyl (photoi), adx=0.25",\
     mypath."/3d_pos_ph_adx1.0_log.txt" u ($2*1e9):7 dt 2 title "3d (photoi), adx=1.0"

set output "cyl_pos_other.pdf"
set title "Cylindrical (no photoi): other effects"
plot mypath."/cyl_pos_adx1.0_log.txt" u ($2*1e9):7 title "adx=1.0",\
     mypath."/cyl_pos_adx1.0_keep_ref_log.txt" u ($2*1e9):7 title "adx=1.0, keep ref.",\
     mypath."/cyl_pos_adx1.0_4vcycle_log.txt" u ($2*1e9):7 title "adx=1.0, 4 v-cycle",\
     mypath."/cyl_pos_adx1.0_halfdt_log.txt" u ($2*1e9):7 title "adx=1.0, half dt",\
     mypath."/cyl_pos_adx1.0_prolong_linear_log.txt" u ($2*1e9):7 title "adx=1.0, prol. linear"

set output "cyl_pos_ph_grid.pdf"
set title "Cylindrical (with photoi.): effect of grid"
plot mypath."/cyl_pos_ph_adx2.0_log.txt" u ($2*1e9):7 title "adx=2.0",\
     mypath."/cyl_pos_ph_adx1.0_log.txt" u ($2*1e9):7 title "adx=1.0",\
     mypath."/cyl_pos_ph_adx0.5_log.txt" u ($2*1e9):7 title "adx=0.5",\
     mypath."/cyl_pos_ph_adx0.25_log.txt" u ($2*1e9):7 title "adx=0.25",\
     mypath."/cyl_pos_ph_unif_6.0_log.txt" u ($2*1e9):7 dt 2 title "unif. 6.0mu",\
     mypath."/cyl_pos_ph_unif_3.0_log.txt" u ($2*1e9):7 dt 2 title "unif. 3.0mu",\
     mypath."/cyl_pos_ph_unif_1.5_log.txt" u ($2*1e9):7 dt 2 title "unif. 1.5mu"

set output "cyl_pos_ph_other.pdf"
set title "Cylindrical (with photoi.): other effects"
plot mypath."/cyl_pos_ph_adx1.0_log.txt" u ($2*1e9):7 title "adx=1.0",\
     mypath."/cyl_pos_ph_adx1.0_keep_ref_log.txt" u ($2*1e9):7 title "adx=1.0, keep ref.",\
     mypath."/cyl_pos_ph_adx1.0_4vcycle_log.txt" u ($2*1e9):7 title "adx=1.0, 4 v-cycle",\
     mypath."/cyl_pos_ph_adx1.0_halfdt_log.txt" u ($2*1e9):7 title "adx=1.0, half dt",\
     mypath."/cyl_pos_ph_adx1.0_prolong_linear_log.txt" u ($2*1e9):7 title "adx=1.0, prol. linear",\
     mypath."/cyl_pos_ph_adx1.0_1step_log.txt" u ($2*1e9):7 title "adx=1.0, 1-step",\
     mypath."/cyl_pos_ph_adx1.0_4step_log.txt" u ($2*1e9):7 title "adx=1.0, 4-step",\
     mypath."/cyl_pos_ph_adx1.0_8step_log.txt" u ($2*1e9):7 title "adx=1.0, 8-step"


