fd = './'
fx = fd.'sb.example30.kd25.0_kb0.0_pu1e-03_rt2.0.R0.e0.'
gx = fd.'sb.example30.kd25.0_kb0.0_pu1e-03_rt2.0.R0.e1.'
hx = fd.'sb.example30.kd25.0_kb0.0_pu1e-03_rt2.0.R0.eX.'

set style line 9 lt 1 lc rgb 'cyan'
array cs[4] = [7, 9, 6, 5]
array lbs[4] = ['M1', 'M2', 'P1', 'P2']
array prs[2] = ['M1-P1', 'M2-P2']
mlw = 1.5
set tics front

set term wxt 1
set view equal xyz
set xyplane 0
set key righ top samplen 0.3
set title 'stage 1 (w/o pairing)' offset 0, -0.5
splot for [cdx=1:4] fx.'ckp' u ($6==(cdx-1) ? $2 : NaN):3:4 w l ls cs[cdx] lw mlw t lbs[cdx]

set term wxt 2
set view equal xyz
set xyplane 0
set key righ top samplen 0.3
set title 'stage 2 (w pairing)' offset 0, -0.5
splot for [cdx=1:4] gx.'ckp' u ($6==(cdx-1) ? $2 : NaN):3:4 w l ls cs[cdx] lw mlw t lbs[cdx]

set term wxt 3
set multiplot layout 1, 2
set xlabel 'Time [10^{3}  frame]'
set xrange [0:]
set xtics 0,0.5,2
set ylabel 'R_{ee} [a]'
set yrange [0:40]
set ytics 0,10,40
set key samplen 0.3 top right
set title 'stage 1 (w/o pairing)' offset 0, -0.5
plot for [cdx=1:4] fx.'log' u ($0/1000):(column(cdx+1)) w l ls cs[cdx] lw mlw t lbs[cdx]
set title 'stage 2 (w pairing)' offset 0, -0.5
plot for [cdx=1:4] gx.'log' u ($0/1000):(column(cdx+1)) w l ls cs[cdx] lw mlw t lbs[cdx]
unset multiplot

set term wxt 4
set multiplot layout 1, 2
set xlabel 'Time [10^{3}  frame]'
set xrange [0:]
set xtics 0,0.5,2
set ylabel 'f_{pair}'
set yrange [-0.02:1.02]
set ytics 0,0.2,1
set key samplen 0.3 bottom right
set title 'stage 1 (w/o pairing)' offset 0, -0.5
plot for [pdx=1:2] fx.'log' u ($0/1000):(column(pdx+5)) w l ls pdx lw mlw t prs[pdx]
set title 'stage 2 (w pairing)' offset 0, -0.5
plot for [pdx=1:2] gx.'log' u ($0/1000):(column(pdx+5)) w l ls pdx lw mlw t prs[pdx]
unset multiplot

set term wxt 5
N = 360
set xrange [-0.5:N-0.5]
set yrange [N-0.5:-0.5]
set size square
set tics front
unset xlabel
unset ylabel
set xtics 45,90,315
set ytics 45,90,315
do for [cdx=1:4] {
  set xtics add(lbs[cdx] cdx*90-45)
  set ytics add(lbs[cdx] cdx*90-45)
}
adx = 1
do for [cdx=1:3] {
  set arrow adx from graph 0, first cdx*90 to graph 1, first cdx*90 ls 8 front nohead
  adx = adx+1
  set arrow adx from first cdx*90, graph 0 to first cdx*90, graph 1 ls 8 front nohead
  adx = adx+1
}
load '../codes/reds.pal'
set logscale cb 10
set cbrange [0.005:1]
set title 'stage 2 (w pairing)' offset 0, -0.5
plot hx.'rc1.5_fq0.px' matrix u 1:2:3 w image notitle




