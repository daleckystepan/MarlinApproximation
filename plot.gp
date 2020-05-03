set terminal qt font "Open Sans"


set multiplot layout 2,2 rowsfirst

#set xrange [-1:0.95]
set yrange [0:5]
plot "plot.dat" using 1:2 with lines title "1/acos", \
             "" using 1:3 with lines title "1/jd marlin", \
             "" using 1:5 with lines title "1/approx", \
             "" using 1:7 with lines title "acosInv", \
             "" using 1:9 with lines title "acosInvQ"

set yrange [-0.1:0.1]
plot 0, \
     "plot.dat" using 1:4 with lines title "1/acos(x) - 1/jd marlin", \
             "" using 1:6 with lines title "1/acos(x) - 1/approx",\
             "" using 1:8 with lines title "1/acos(x) - acosInv", \
             "" using 1:10 with lines title "1/acos(x) - acosInvQ"
              

unset yrange
set xrange [0.995:]
set yrange [10:90]
plot "plot.dat" using 1:2 with lines title "1/acos", \
             "" using 1:3 with lines title "1/jd marlin", \
             "" using 1:5 with lines title "1/approx", \
             "" using 1:7 with lines title "acosInv", \
             "" using 1:9 with lines title "acosInvQ"

set xrange [0.999:]
set yrange [-5:5]
plot 0, \
     "plot.dat" using 1:4 with lines title "1/acos(x) - 1/jd marlin", \
             "" using 1:6 with lines title "1/acos(x) - 1/approx",\
             "" using 1:8 with lines title "1/acos(x) - acosInv", \
             "" using 1:10 with lines title "1/acos(x) - acosInvQ"

unset multiplot