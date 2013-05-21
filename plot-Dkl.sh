#!/bin/bash

make_plts () {
  gnuplot -e "set terminal dumb;
    set title \"dD\/dmx\";
    plot \"$1\" using (\$8-\$8*\$10) w dot"

  gnuplot -e "set terminal dumb;
    set title \"dD\/dmy\";
    plot \"$1\" using (\$11-\$12*\$13) w dot"

  gnuplot -e "set terminal dumb;
    set title \"dD\/dBxx\";
    plot \"$1\" using (\$14-\$15*\$16) w dot"

  gnuplot -e "set terminal dumb;
    set title \"dD\/dByy\";
    plot \"$1\" using (\$17-\$18*\$19) w dot"
}

make_plts $1 | less
