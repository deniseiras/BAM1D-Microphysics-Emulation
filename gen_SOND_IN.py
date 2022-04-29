import random

for year in [ '2014', '2015']:

    for expnum in range(1, 11):
        f = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_{}/SOND_IN_original'.format(year))
        fw = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_{}/SOND_IN_NEW_{}'.format(year, expnum), 'w')

        for l in f:
            line = l
            spl = line.split()
            linew = '{:>13.6f}'.format(float(spl[0]))
            for i in range(1, len(spl)):
                linew += '{:>13.6f}'.format(float(spl[i]) + float(spl[i]) * random.uniform(-0.05, 0.05))
            linew += '\n'
            fw.write(linew)

        f.close()
        fw.close()

