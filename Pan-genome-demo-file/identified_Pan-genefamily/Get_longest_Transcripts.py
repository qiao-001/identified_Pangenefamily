####仅适用于.区分的转录本
####仅适用于.区分的转录本
from pathlib import Path

p = Path("D:/13.pan_genome/5.Camelina_sativa_LBDpan/LBD/4.3OrthoFinder/13protein_")
def select_longest(seq_fa,longest_fa):
    #初始化空字典用于存放序列，字典的key是基因名，value存放转录本和序列
    #格式为>{gene1:[>gene1.1\nabd\n>gene1.2\neee\nerr\n....]}
    seq_dict = {}
    with open(seq_fa,'rt') as f1:
        for eachline in f1:  
            if eachline[0] == '>':
                gene_name = eachline[1:].split('.')[0]#更改此处.可以选其他符号区分的转录本
                if gene_name not in seq_dict.keys():
                    seq_dict[gene_name] = eachline #此处不加strip(),方便用\n区别转录本和序列
                else:
                    seq_dict[gene_name] += eachline
            elif len(eachline.strip()) > 0:
                seq_dict[gene_name] += eachline
            

    with open(longest_fa,'wt') as f2:
        for i in seq_dict.keys():
            transcripts_list = seq_dict[i].split('>')[1:] #此处用切片去除列表的第一个空字符,得到列表 ['gene1\n12\n3456\n', 'gene2\n123\n45\n', 'gene3\n1234\n56789\n']
            seq_list= []
            for each in transcripts_list:
                each = each.split('\n')#单独提取序列
                seq = ''.join(each[1:])
                seq_list.append(seq)

            length = list(map(len,seq_list))
            index = length.index(max(length))
            f2.write(f'>{transcripts_list[index]}')

def search_fa_longest():
    for file in p.iterdir():
        if file.suffix == '.fa':
            outfile = file.parent / (file.stem + '_L.fa')
            select_longest(file,outfile)

search_fa_longest()      