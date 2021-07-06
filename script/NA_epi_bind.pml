set ray_shadows, 0
bg_color white
set valence, 0
fetch 2AEP, N2_head
fetch 4GZO, N2


create NA, chain A and N2_head

select epi, resi 328-329+344+367-370 and NA

hide all
color white, NA
show cartoon,NA

show sticks, resi 328+344+369 and (not name c+n+o)
color atomic, epi
util.cbaw

distance epi_ds, /NA/A/A/GLU`344/OE1, /NA/A/A/LYS`369/NZ
distance epi_ds, /NA/A/A/GLU`344/OE2, /NA/A/A/LYS`328/NZ
distance epi_ds, /NA/A/A/LYS`328/NZ, /NA/A/A/LYS`369/NZ
color black, epi_ds
hide label
set_view (\
    -0.446687430,    0.555872440,   -0.701050460,\
     0.885415733,    0.162127033,   -0.435602307,\
    -0.128482103,   -0.815299749,   -0.564598978,\
    -0.000258662,   -0.000294938,  -45.616333008,\
   117.913421631,   55.264221191,   35.213424683,\
    32.371498108,   58.991981506,  -20.000000000 )
ray
png ../graph/NA_epi_bind1.png
hide all
select epi2, resi 328-329+344+367-370 and N2
color white, N2
show cartoon,N2

show sticks, resi 328+344 and (not name c+n+o)
color atomic, epi2
util.cbaw
distance epi_ds3, /N2/A/A/GLU`344/OE2, /N2/A/A/LYS`328/NZ
color black, epi_ds3
hide label

set_view (\
    -0.980715215,    0.017208850,    0.194672048,\
    -0.165143222,    0.459670812,   -0.872595012,\
    -0.104503043,   -0.887918770,   -0.447963715,\
     0.000008106,    0.000049744,  -53.938377380,\
    -5.118977547,   11.772624016,   11.733413696,\
    44.771537781,   63.123916626,  -20.000000000 )
ray
png ../graph/NA_epi_bind2.png