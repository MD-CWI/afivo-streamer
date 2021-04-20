mkdir output/mesh/

./streamer standard.cfg -output%name="output/mesh/adx02" -refine_adx=0.2

./streamer standard.cfg -output%name="output/mesh/adx1" -refine_adx=1.0

./streamer standard.cfg -output%name="output/mesh/adx2" -refine_adx=2.0

#./streamer standard.cfg -output%name="output/adx01" -refine_adx=0.1

#./streamer standard.cfg -output%name="output/adx05_fac1.2" -refine_adx=0.5 -refine_adx_fac=1.2

#./streamer standard.cfg -output%name="output/adx05_defe5" -refine_adx=0.5 -derefine_dx=1e-5
