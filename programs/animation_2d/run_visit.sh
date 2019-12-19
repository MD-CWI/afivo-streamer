#!/bin/bash

# Author: Jannis Teunissen
# Instructions: uncomment the steps that you want to perform below

mkdir -p output
mkdir -p pngs

# Run the simulation
# ./streamer cyl.cfg

# Produce figures with visit
visit -nowin -cli -s visit_plot_2d.py "output/cyl_*.silo" electric_fld \
      -outdir pngs -cmin 0 -cmax 1.5e7 -nodal -bbox 0. 2e-3 0. 16e-3 \
      -contour phi -contour_nlevels 17 -contour_cmin 0 -contour_cmax 3.2e4 \
      -contour_linewidth 2 -contour_rgb 255 255 255 &

visit -nowin -cli -s visit_plot_2d.py "output/cyl_*.silo" rhs \
      -outdir pngs -cmin " -2e11" -cmax 2e11 -nodal -bbox 0. 2e-3 0. 16e-3 \
      -ct difference &

visit -nowin -cli -s visit_plot_2d.py "output/cyl_*.silo" e \
      -outdir pngs -cmin 0 -cmax 2e20 -nodal -bbox 0. 2e-3 0. 16e-3 &

visit -nowin -cli -s visit_plot_2d.py "output/cyl_*.silo" N2_C3 \
      -outdir pngs -cmin 0 -cmax 1e20 -nodal -bbox 0. 2e-3 0. 16e-3 \
      -ct inferno &

visit -nowin -cli -s visit_plot_2d.py "output/cyl_*.silo" N2_plus \
      -outdir pngs -cmin 1e40 -cmax 2e40 -nodal -bbox 0. 2e-3 0. 16e-3 \
      -ct xray -mesh blk_mesh &

wait

# Trim all the images (here done with GNU parallel)
# parallel convert -trim {} pngs/small_{/} ::: pngs/cyl_*.png

# Combine the images
# for i in {0..220}; do
#     convert \
#         \( "pngs/small_cyl_electric_fld_$i.png" -flop \) "pngs/small_cyl_electric_fld_$i.png" \
#         \( "pngs/small_cyl_e_$i.png" -flop \) "pngs/small_cyl_e_$i.png" \
#         \( "pngs/small_cyl_rhs_$i.png" -flop \) "pngs/small_cyl_rhs_$i.png" \
#         \( "pngs/small_cyl_N2_C3_$i.png" -flop \) "pngs/small_cyl_N2_C3_$i.png" \
#         +append "pngs/combined_$i.png"

#     convert \
#         \( "pngs/small_cyl_electric_fld_$i.png" -flop \) "pngs/small_cyl_electric_fld_$i.png" \
#         \( "pngs/small_cyl_e_$i.png" -flop \) "pngs/small_cyl_e_$i.png" \
#         \( "pngs/small_cyl_rhs_$i.png" -flop \) "pngs/small_cyl_rhs_$i.png" \
#         \( "pngs/small_cyl_N2_C3_$i.png" -flop \) "pngs/small_cyl_N2_C3_$i.png" \
#         \( "pngs/small_cyl_N2_plus_$i.png" -flop \) "pngs/small_cyl_N2_plus_$i.png" \
#         +append "pngs/combined_mesh_$i.png"
# done

# Convert them to clips
# ffmpeg -y -i pngs/combined_%d.png -c:v libx264 -crf 22 -pix_fmt yuv420p -vf "scale=1600:-1" cyl_combined.mp4
# ffmpeg -y -i pngs/combined_%d.png -c:v libx264 -crf 22 -pix_fmt yuv420p-vf "scale=800:-1" cyl_combined_small.mp4
# ffmpeg -y -i pngs/combined_mesh_%d.png -c:v libx264 -crf 22 -pix_fmt yuv420p -vf "scale=1600:-1" cyl_combined_mesh.mp4
# ffmpeg -y -i pngs/combined_mesh_%d.png -c:v libx264 -crf 22 -pix_fmt yuv420p -vf "scale=800:-1" cyl_combined_mesh_small.mp4
# ffmpeg -y -i pngs/combined_%d.png -c:v libvpx -b:v 1M -pix_fmt yuv420p -vf "scale=1920:-1" cyl_combined_small.webm
