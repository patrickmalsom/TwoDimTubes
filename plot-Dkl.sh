#!/bin/bash

# ======== questions for Lin Clon ========
# how to deal with a null for pltnum
# how to querry for all 2WellTube* files and allow user to choose one


echo "This script generates terminal plots of configurations for"
echo "the Tubes HMC routine for a 2D potential"

pltnum=1
while [ $pltnum -gt 0 ]; 
do 
  echo " "
  echo "========================================================================="
  echo "Enter 0 to exit"
  echo "========================================================================="
  echo "Available Plots: "
  echo " (1): input x path     (inputPos.dat)"
  echo " (2): input y path     (inputPos.dat)"
  echo " (3): input mx         (inputMean.dat)"
  echo " (4): input dmx/dt     (inputMean.dat)"
  echo " (5): input d^2mx/dt^2 (inputMean.dat)"
  echo " (6): input Bxx        (inputB.dat)"
  echo " (7): input Byy        (inputB.dat)"
  echo " (8): input Bxy        (inputB.dat)"
  echo " (9): output mx        (E[x])"
  echo "(10): output my        (E[y])"
  echo "(11): output dD/dmx    (E[(Linvdgz -z)*I]-E[Linvdgz-z]*E[I])"
  echo "(12): output dD/dmy    (E[(Linvdgz -z)*I]-E[Linvdgz-z]*E[I])"
  echo "(13): output dD/dBxx   (E[zz*i]-E[zz]*E[I])"
  echo "(14): output dD/dByy   (E[zz*i]-E[zz]*E[I])"
  echo "(15): output dD/dBxy   (E[zz*i]-E[zz]*E[I])"
  echo "(16): output G         (E[G])"
  echo "========================================================================="

  echo "Please enter a plot number, followed by [Enter]:"
  #read the user input
  read pltnum
  echo $pltnum

  case $pltnum in

    "") exit 
       ;;
    0) exit 
       ;;

    # plot the input path (termalized)
    1) gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputPosx\";\
       plot \"inputPos.dat\" using (\$1) w line" 
       ;;
    2) gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputPosx\";\
         plot \"inputPos.dat\" using (\$2) w line" 
       ;;

    # plot the input Mean x
    3) gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputMeanx\";\
         plot \"inputMean.dat\" using (\$1) w line"
       ;;
    4) gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputdMeanx\/dt\";\
         plot \"inputMean.dat\" using (\$2) w line"
       ;;
    5) gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputd2Meanx\/dt2\";\
         plot \"inputMean.dat\" using (\$3) w line"
       ;;

    # plot the input B
    6) gnuplot -e "set terminal dumb nofeed 100 40; set title \"input Bxx\";\
         plot \"inputB.dat\" using (\$1) w line"
       ;;
    7) gnuplot -e "set terminal dumb nofeed 100 40; set title \"input Byy\";\
         plot \"inputB.dat\" using (\$2) w line"
       ;;
    8) gnuplot -e "set terminal dumb nofeed 100 40; set title \"input Bxy\";\
         plot \"inputB.dat\" using (\$3) w line"
       ;;

    #plot the mean path (regular average)
    9) gnuplot -e "set terminal dumb nofeed 100 40; set title \"Mean X Path\";\
         plot \"$1\" using (\$3) w line"
       ;;
    10) gnuplot -e "set terminal dumb nofeed 100 40; set title \"Mean Y Path\";\
         plot \"$1\" using (\$4) w line"
       ;;

    #plot the mean kl distance
    11) gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dmx\";\
         plot \"$1\" using (\$8-\$9*\$10) w line"
       ;;
    12) gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dmy\";\
         plot \"$1\" using (\$11-\$12*\$13) w line"
       ;;

    #plot the tube kl distance
    13) gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dBxx\";\
          plot \"$1\" using (-\$14+\$15*\$16) w line"
       ;;
    14) gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dByy\";\
          plot \"$1\" using (-\$17+\$18*\$19) w line"
       ;;
    15) gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dByy\";\
          plot \"$1\" using (-\$20+\$21*\$22) w line"
       ;;

    #Plot G 
    16) gnuplot -e "set terminal dumb nofeed 100 40; set title \"E[G]\";\
          plot \"$1\" using (\$25) w line"
       ;;
  esac

  read -p "Press [Enter] to make another"

done

