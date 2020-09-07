#!/bin/bash
cd kid_a_setup 

export LD_LIBRARY_PATH=..:bin:$LD_LIBRARY_PATH 
export LD_PRELOAD=ptrutil.so 
export PYTHONPATH=/mnt/local/pdziekan/usr/local/lib/python3/dist-packages/

for sd_conc in 100
do
  for sstp_cond in 10 
  do
    for sstp_coal in 100 
    do
      for ntot in 50e6
      #for ntot in 50e6 150e6 300e6
      do
        for spinup in 0
        #for spinup in 0 1e8
        do
          RAIN="NORAIN"
          if [ ${spinup} == "0" ] ; then
            RAIN="RAIN"
          fi
          for wctrl in 2
          #for wctrl in 2 3
          do
            if [ ${wctrl} == "2" ] ; then
              sed -i 's/wctrl(1)=3/wctrl(1)=2/g' namelists/kida_icmw1D_libcloud_lgr.nml
            fi
            if [ ${wctrl} == "3" ] ; then
              sed -i 's/wctrl(1)=2/wctrl(1)=3/g' namelists/kida_icmw1D_libcloud_lgr.nml
            fi
            OUTDIR=/mnt/local/pdziekan/wyniki/kid-a/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl${wctrl}_ntot${ntot}_${RAIN}_HALL && mkdir ${OUTDIR}
            for iter in {1..40}
            do
              mkdir ${OUTDIR}/${iter} 
              FILEOUT=${OUTDIR}/${iter} python3 ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --backend="serial" --sd_conc=$sd_conc --n_tot=$ntot --spinup_rain=$spinup &
              sleep 2
            done
            wait
            for iter in {41..80}
            do
              mkdir ${OUTDIR}/${iter} 
              FILEOUT=${OUTDIR}/${iter} python3 ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --backend="serial" --sd_conc=$sd_conc --n_tot=$ntot --spinup_rain=$spinup &
              sleep 2
            done
            wait
            for iter in {81..120}
            do
              mkdir ${OUTDIR}/${iter} 
              FILEOUT=${OUTDIR}/${iter} python3 ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --backend="serial" --sd_conc=$sd_conc --n_tot=$ntot --spinup_rain=$spinup &
              sleep 2
            done
            wait
          done
        done
      done
    done
  done
done

echo "done with SD100 ensemble"
