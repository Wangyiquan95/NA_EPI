set ray_shadows, 0
bg_color white
set valence, 0
fetch 4gzq, N2_head



select NA, chain A and N2_head
select active_site1, resi 118+119+134+151+152+178+179+224-227+276+277+292+371+406 and NA
select SIA, /N2_head/C/A/SIA
select epi, resi 328-329+344+367-370 and NA

hide all


color white, NA
show cartoon,NA
show sticks, SIA
show spheres, epi and name ca
set sphere_scale, 1
color lightblue, epi
color yellow, SIA
color deepsalmon, active_site
cartoon tube
util.cnc all

set_view (\
    -0.952349842,   -0.097421058,    0.289029270,\
    -0.258540034,    0.760601342,   -0.595517278,\
    -0.161823049,   -0.641867101,   -0.749545634,\
    -0.000022657,   -0.000007644,  -85.927253723,\
     3.029790163,   16.614337921,   10.380401611,\
    71.233283997,  100.616760254,  -20.000000000 )
ray
png ../graph/NA_epi_zoom.png
