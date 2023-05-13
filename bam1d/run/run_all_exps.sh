# Deixar igual ao runbam !!!
# So serve pra diferenciar os diretorios de saída
#

levs='64'
deltat='90'

#ILCON='LSC'
# csv_prefix='none'
#ILCON='HGFS'
# csv_prefix='ferrier'
ILCON='HUMO_NN'   # HUMO_NN para NN !!!!!
csv_prefix='hug_morr_'

ISCON='SHALLOFF'
ICCON='DEEPOFF'
#ISCON='TIED'
#ICCON='GRE'

ISWRAD='CRD'
ILWRAD='RRTMG' 	#(RRTMG HRS HRS) 

# refere-se a de quantos em quantos timesteps se grava o arquivo csv (Micro_HugMorr.f90)
# Ex saveeach='2' - salvando a cada 2 timesteps de 360 (base = 360)
# Isso faz com que, ao variar deltat, a quantidade de dados salva seja fixa para um mesmo numero de niveis
saveeach='2'

# gen data
# for expy in "2014" "2015"; do

# test NN para treinamento 2014 e parte de 2015. Testes, 2015
for expy in "2014" "2015"; do

  direxp=/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_${expy}

  rm ${direxp}/SOND_IN
  ln -s ${direxp}/SOND_IN_original ${direxp}/SOND_IN
  ./run_expname_expdir.sh IOP${expy}_${ISWRAD}${ILWRAD}_${ILCON}__${ISCON}_${ICCON}_${deltat}_savedeach${saveeach}_${levs} GOAMAZON_${expy}

done

# # Data Augmentation
# #
# for expn in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"; do 
#   echo ${direxp}/SOND_IN_NEW_${expn}
#   rm ${direxp}/SOND_IN
#   ln -s ${direxp}/SOND_IN_NEW_${expn} ${direxp}/SOND_IN
#   ./run_expname_expdir.sh IOP${expy}_CRDRRTMG_${ILCON}_${deltat}${levs}_SOND_NEW_${expn} GOAMAZON_${expy}
# done


