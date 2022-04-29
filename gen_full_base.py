#
# Um ano tem 810880−(80*7*28) = 795200 linhas
# corresponde a 795200÷28 = 28400 registros
# aproximadamente 200 MB
#
levels=18
spin_times = 80*7*levels  # saving 80 times / day => discard 7 days from start
max_ensembles=5
path='/media/denis/dados/_COM_BACKUP/NN_BAM1D/bam1d_data'

for inp_out in ['in', 'out']:
    fw = open('{}/IOP_2014-2015_CRDRRTMG_360{}_FULL{}x___hug_morr_{}puts.csv'.format(path, levels, max_ensembles, inp_out), 'w')
    first = True
    for year in [ '2014', '2015']:
        for expnum in range(1, max_ensembles+1):
            path2 = '{}/IOP{}_CRDRRTMG_360{}_SOND_NEW_{}'.format(path, year, levels, expnum)
            f = open('{}/hug_morr_{}puts.csv'.format(path2, inp_out))

            header = f.readline()

            if first:
                fw.write(header)
                first = False

            # skip spinup lines
            for spinup in range(spin_times):
                _ = f.readline()

            for line in f:
                fw.write(line)
            f.close()

    fw.close()


