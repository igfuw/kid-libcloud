#!/bin/bash
cd kid_a_setup 


export LD_LIBRARY_PATH=..:bin:$LD_LIBRARY_PATH 
export LD_PRELOAD=ptrutil.so 
export PYTHONPATH=/mnt/local/pdziekan/lib64/python2.7/site-packages/

for sd_conc in 40
do
  for sstp_cond in 10 
  do
    for sstp_coal in 10 
    do
      for ntot in 50e6 150e6 300e6
      do
        for spinup in 0 1e8
        do
          RAIN="NORAIN"
          if [ ${spinup} == "0" ] ; then
            RAIN="RAIN"
          fi
          for wctrl in 2 3
          do
            if [ ${wctrl} == "2" ] ; then
              sed -i 's/wctrl(1)=3/wctrl(1)=2/g' namelists/kida_icmw1D_libcloud_lgr.nml
            fi
            if [ ${wctrl} == "3" ] ; then
              sed -i 's/wctrl(1)=2/wctrl(1)=3/g' namelists/kida_icmw1D_libcloud_lgr.nml
            fi
            OUTDIR=/mnt/local/pdziekan/wyniki/kid-a/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl${wctrl}_ntot${ntot}_${RAIN} && mkdir ${OUTDIR}
            for iter in {1..2}
            do
              mkdir ${OUTDIR}/${iter} &
              FILEOUT=${OUTDIR}/${iter} python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --backend="serial" --sd_conc=$sd_conc --n_tot=$ntot --spinup_rain=$spinup &
            done
            wait
          done
        done
      done
    done
  done
done

echo "done with SD40 ensemble"

#mkdir output
#for sd_conc in 10000 #100 1000 10000
#do
#  for sstp_cond in 10 100
#  do
#    for sstp_coal in 10 100
#    do
#      sed -i 's/wctrl(1)=2/wctrl(1)=3/g' namelists/kida_icmw1D_libcloud_lgr.nml
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=5e7 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot50_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=1.5e8 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot150_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=3e8 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot300_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=5e7 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot50 && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=1.5e8 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot150 && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=3e8 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl3_ntot300 && mkdir output
#      sed -i 's/wctrl(1)=3/wctrl(1)=2/g' namelists/kida_icmw1D_libcloud_lgr.nml
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=5e7 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot50_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=1.5e8 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot150_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=3e8 --spinup_rain=1e6
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot300_NoCoalNoSedi && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=5e7 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot50 && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=1.5e8 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot150 && mkdir output
#      env python ../kid_1D.py --sstp_coal=$sstp_coal --sstp_cond=$sstp_cond --sd_conc=$sd_conc --n_tot=3e8 --spinup_rain=0
#      mv output $HOME/praca/wyniki/KiD-A/1D/output_SdConc${sd_conc}_SstpCond${sstp_cond}_SstpCoal${sstp_coal}_wctrl2_ntot300 && mkdir output
#    done
#  done
#done

