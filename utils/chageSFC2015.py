
f = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_2014-2015/SFC_ASCII_2015')
fw = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_2014-2015/SFC_ASCII_2015_for_2014', 'w')



for i, l in enumerate(f):
    if i==1: continue
    line = l
    spl = line.split(maxsplit=1)
    spl[0] = str('{:.6f}'.format(float(spl[0]) + 364.861115 + (365/13136)))
    line = '   ' + spl[0] + '   ' + spl[1]
    print(line)
    fw.write(line)

