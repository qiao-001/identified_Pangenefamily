from pathlib import Path

files = Path("D:/13.pan_genome/5.Camelina_sativa_LBDpan/LBD/4.2OrthoFinder")#change this file path

def id_simplify(i,outfile):
    with open(i,'rt') as f1:
        with open(outfile,'wt') as f2:
            for eachline in f1:
                if eachline[0] == '>':
                    f2.write(eachline.strip().split()[0])
                    f2.write('\n')
                else:
                    f2.write(eachline)
                    
def search_fa_simplify():
    for i in files.iterdir():
        if i.suffix == '.fa':
            outfile = i.parent / (i.stem + '_sim.fa')
            id_simplify(i,outfile)

search_fa_simplify()