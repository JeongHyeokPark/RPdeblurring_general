#!/bin/bash


iter=0

isbatch=1

if [ $isbatch -eq 1 ]; then
  while :
  do
  chknum=0
  cd /home/jhpark/work/deblurring_spirit/merging/results;

  if [ $iter -eq 0 ]; then
    rm /home/jhpark/public_html/files/deblurring_spirit_eff/new_eff/*png;
    rm /home/jhpark/work/deblurring_spirit/merging/results/tmp_files/*.root;
    rm /home/jhpark/work/deblurring_spirit/merging/results/*.root;
    rm /home/jhpark/work/deblurring_spirit/merging/*.root;

    cd /home/jhpark/work/deblurring_spirit/merging/;
    ./data_production.exe Input.txt;
    cp output_dist_Sn108.root input.root;
    ./submit_batch_jobs.sh $iter;
    cd /home/jhpark/work/deblurring_spirit/merging/results;
  fi

  # Check if there are correct numbers of the output files in the results directory
  while :
    do
    number=$(ls output_* | wc -l)
    if [ $(($chknum % 100)) -eq 0 ]; then
    echo "iteration $iter checking the output numbers $number"
    fi

    if [ $number -eq 1000 ]; then
      ((iter++))
    break
    fi

    ((chknum++))
    sleep 1
  done

  # If the condor jobs are finished, start RL algorithm
  echo "start $iter"

  # Copy the files created from the script, accumulate_and...
  ./make_and_copy.sh $iter;

  # Start RL algorithm
  cd /home/jhpark/work/deblurring_spirit/merging/;
  ./do_deblurring.exe Input.txt ReadyToDeblur_$iter.root;
  root -q -b draw.C\($iter\);

  if [ $iter -eq 20 ]; then
  #if [ $iter -eq 50 ]; then
    break
  fi

  # Submit the condor running files
  ./submit_batch_jobs.sh $iter;

  done

  ./error_calculation.exe Input.txt input.root
fi





if [ $isbatch -eq 0 ]; then
  while :
  do

  cd /home/jhpark/work/deblurring_spirit/merging;

  if [ $iter -eq 0 ]; then
    rm /home/jhpark/public_html/files/deblurring_spirit_eff/new_eff/*png;
    rm /home/jhpark/work/deblurring_spirit/merging/results/tmp_files/*.root;
    rm /home/jhpark/work/deblurring_spirit/merging/*.root;
    ./data_production.exe Input.txt;
    cp output_dist_Sn108.root input.root;
  fi

  ((iter++))
  ./call_and_save.exe Input.txt input.root $iter;

  echo "start $iter"

  # Start RL algorithm
  ./do_deblurring.exe Input.txt output_$iter.root;
  root -q -b draw.C\($iter\);

  if [ $iter -eq 50 ]; then
    break
  fi

  done

  ./error_calculation.exe Input.txt input.root
fi
