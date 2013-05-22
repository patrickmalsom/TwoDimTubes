#!/bin/bash

make_plts () {
#plot the input path
gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputPosx\"; plot \"inputPos.dat\" using (\$1) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputPosy\"; plot \"inputPos.dat\" using (\$2) w line"

#plot the input mean x path
gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputMeanx\"; plot \"inputMean.dat\" using (\$1) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputdMeanx\/dt\"; plot \"inputMean.dat\" using (\$2) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"inputd2Meanx\/dt2\"; plot \"inputMean.dat\" using (\$3) w line"

#plot the mean path (regular average)
gnuplot -e "set terminal dumb nofeed 100 40; set title \"Mean X Path\"; plot \"$1\" using (\$3) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"Mean Y Path\"; plot \"$1\" using (\$4) w line"

#plot the mean kl distance
gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dmx\"; plot \"$1\" using (\$8-\$9*\$10) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dmy\"; plot \"$1\" using (\$11-\$12*\$13) w line"

#plot the tube kl distance
gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dBxx\"; plot \"$1\" using (-\$14+\$15*\$16) w line"
gnuplot -e "set terminal dumb nofeed 100 40; set title \"dD\/dByy\"; plot \"$1\" using (-\$17+\$18*\$19) w line"

}

make_plts $1 | less
