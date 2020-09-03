#!/bin/bash
cd kid_a_setup 


export LD_LIBRARY_PATH=..:bin:$LD_LIBRARY_PATH 
export LD_PRELOAD=ptrutil.so 
export PYTHONPATH=/mnt/local/pdziekan/usr/local/lib/python3/dist-packages/

gpu_number=4 # number of gpus to be used by this simulation

for run_no in {1..10}
do
  for spinup in 0 1e8
  do
    RAIN="NORAIN"
    if [ ${spinup} == "0" ] ; then
      RAIN="RAIN"
    fi
    for sstp_cond in 10
    do
      for sstp_coal in 100
      do
        for sd_conc in 10000
        #for sd_conc in 10000 0
        do
          for ntot in 50e6 150e6 300e6
          do
            for wctrl in 2 3
            do
              if [ ${wctrl} == "2" ] ; then
                sed -i 's/wctrl(1)=3/wctrl(1)=2/g' namelists/kida_icmw1D_libcloud_lgr.nml
              fi
              if [ ${wctrl} == "3" ] ; then
                sed -i 's/wctrl(1)=2/wctrl(1)=3/g' namelists/kida_icmw1D_libcloud_lgr.nml
              fi

              sd_const_multi=0
              #if [ ${sd_conc} == "0" ] ; then
              #  if [ ${ntot} == "50e6" ] ; then
              #    sd_const_multi=125000
              #  fi
              #  if [ ${ntot} == "150e6" ] ; then
              #    sd_const_multi=375000
              #  fi
              #  if [ ${ntot} == "300e6" ] ; then
              #    sd_const_multi=750000
              #  fi
              #fi

              OUTDIR=/mnt/local/pdziekan/wyniki/kid-a/1D/output_SdConc${sd_conc}_SdConstMulti${sd_const_multi}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl${wctrl}_ntot${ntot}_${RAIN} && mkdir ${OUTDIR}

              # look for an idle gpu and run on it
              iter=-1 # simulation iterator
              while :
              do
                iter=$(($iter+1))
                gpu_id=$(($iter%$gpu_number))
                if [[ $(nvidia-smi -i ${gpu_id}  --query-compute-apps=pid --format=csv,noheader) ]]; then
                  :
                else
                  echo "gpu ${gpu_id} is idle, running on it"
                  #CUDA_DEVICE_ORDER set to PCI_BUS_ID to get consistency of CUDA_VISIBLE_DEVICES setting with nvidia-smi output
                  CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=${gpu_id} FILEOUT=${OUTDIR}/${run_no} python3 ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --backend="CUDA" --sd_conc=$sd_conc --sd_const_multi=$sd_const_multi --n_tot=$ntot --spinup_rain=$spinup &
                  sleep 10 # startup make sime time
                  break
                fi
                sleep 2
              done
            done
          done
        done
      done
    done
  done
done

echo "done"
