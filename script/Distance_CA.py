
from pymol import cmd

N2=cmd.load('script/2aep.cif','N2')
resi_ls=['/N2/A/A/LYS`328/CA','/N2/A/A/ASN`329/CA','/N2/A/A/GLU`344/CA','/N2/A/A/SER`367/CA','/N2/A/A/GLU`368/CA','/N2/A/A/LYS`369/CA','/N2/A/A/PHE`370/CA']
output='result/CA_distance.tsv'


def Distance_CA(resi_ls,output):
    outfile = open(output, 'w')
    outfile.write('pair' + "\t" + 'distance' + "\n")
    # compute distance between CA atoms on antigenic region
    for resi1 in resi_ls:
        for resi2 in resi_ls:
            dst=cmd.get_distance(resi1,resi2)
            pair1=resi1.rsplit('/')[4]
            pair2=resi2.rsplit('/')[4]
            pair=pair1+'-'+pair2
            outfile.write(str(pair) + "\t" + str(dst) + "\n")
    outfile.close()

Distance_CA(resi_ls,output)


