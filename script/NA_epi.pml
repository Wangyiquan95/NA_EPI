set ray_shadows, 0
bg_color white

fetch 4gzq, N2_head

#Smooth surface quick (blob)
set solvent_radius, 2
set valence, 0
symexp tetramer, N2_head, (N2_head),2.5
select NA, chain A and N2_head
select active_site, resi 118+119+134+151+152+178+179+224-227+276+277+292+371+406 and NA
select SIA, /N2_head/C/A/SIA
select epi, resi 328-329+344+367-370 and NA


color white, NA
color lightblue, epi
color yellow, SIA
color deepsalmon, active_site

color grey60, tetramer01 or tetramer03 or tetramer02


hide all

show surface, epi
show surface, active_site
show sticks, SIA
show surface, NA
show surface, tetramer01 or tetramer03 or tetramer02
set_view (\
    -0.944328487,   -0.223620683,    0.241316736,\
    -0.322693139,    0.486636400,   -0.811821342,\
     0.064106733,   -0.844496369,   -0.531706274,\
     0.000003994,    0.000030458, -266.506652832,\
     2.673481941,   38.424980164,   14.691686630,\
   214.573181152,  318.440856934,  -20.000000000 )
 ray

 png ../graph/NA_epi.png