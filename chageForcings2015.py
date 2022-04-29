
f = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_2014-2015/FORCINGS_ASCII 2015')
fw = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_2014-2015/FORCINGS_ASCII 2015_for_2014', 'w')

for l in f:
    line = l
    spl = line.split()
    if len(spl) == 5:
        spl = line.split(maxsplit=2)
        spl[0] = str(int(spl[0]) + 26271)
        line = spl[0] + ' ' + spl[1]
        print(line)
    fw.write(line)

