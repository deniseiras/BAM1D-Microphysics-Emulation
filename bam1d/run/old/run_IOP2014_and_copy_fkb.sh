
# copia modelo recem baixado e convertido no jupyter local (seção FKB Test)
cp /home/denis/_COM_BACKUP/NN_BAM1d/bam1d_physics_nn/fkb_model.txt /home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_2014/

now=$(date '+%Y-%m-%d_%H_%M_%S')
./run_bam GOAMAZON_2014 &> out_IOP2014_${now}.txt