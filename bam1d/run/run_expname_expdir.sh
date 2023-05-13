now=$(date '+%Y-%m-%d_%H_%M_%S')
expname=$1
modelin_dirname=$2
csv_prefix=$3

dir_csv=/media/denis/dados/_COM_BACKUP/NN_BAM1D/bam1d_data/${expname}

if [ -d "$dir_csv" ]; then
  echo "Dir exists: ${dir_csv} . exiting ..."
	exit -1
fi

fileout=out_${expname}_${now}.txt
touch ${fileout}
# tail -f ${fileout} &
./run_bam ${modelin_dirname} &> /dev/null # ${fileout} 

# ./run_bam ${modelin_dirname}

rm -rf ../model/dataout/model_dataout/${expname}
mv ../model/dataout/Exp_01 ../model/dataout/model_dataout/${expname}



rm -rf ${dir_csv}
mkdir ${dir_csv}

mv /home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/${csv_prefix}*puts.csv ${dir_csv}
