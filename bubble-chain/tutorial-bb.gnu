fd = './'
fx = fd.'bb.example30.s0_kd25.0_dp0.7937.R0.e0.'
gx = fd.'bb.example30.s1_kd25.0_dp0.7937.R0.e0.'
ux = fd.'bb.example30.s0_kd25.0_dp0.7937.R0.eX.'
vx = fd.'bb.example30.s1_kd25.0_dp0.7937.R0.eX.'

set style line 9 lt 1 lc rgb 'cyan'
array cs[4] = [7, 9, 6, 4]
array lbs[4] = ['M1', 'M2', 'P1', 'P2']
array prs[2] = ['M1-P1', 'M2-P2']
array dts[4] = [1,1,3,3]
mlw = 1.5
mps = 0.7
set tics front

set term wxt 1
set view equal xyz
set xyplane 0
set key righ top samplen 0.3
splot for [cdx=1:3:2] fx.'ckp' u ($6==(cdx-1)/2 ? $2 : NaN):3:4 w l ls cs[cdx] lw mlw t lbs[cdx], \
      fx.'ckp' u ($6==0 & $5==1 ? $2 : NaN):3:4 w p ls 8 lw mlw pt 7 ps mps t 'button'

set term wxt 2
set view equal xyz
set xyplane 0
set key righ top samplen 0.3
splot for [cdx=2:4:2] gx.'ckp' u ($6==(cdx-2)/2 ? $2 : NaN):3:4 w l ls cs[cdx] lw mlw t lbs[cdx], \
      gx.'ckp' u ($6==0 & $5==1 ? $2 : NaN):3:4 w p ls 8 lw mlw pt 7 ps mps t 'button'

set term wxt 3
set multiplot layout 1, 2
set xlabel 'Time [10^{3}  frame]'
set xrange [0:]
set xtics 0,0.5,2
set ylabel 'R_{ee} [a]'
set yrange [0:40]
set ytics 0,10,40
set key samplen 0.3 top right
set title 'M1 - P1' offset 0, -0.5
plot for [cdx=1:3:2] fx.'log' u ($0/1000):(column((cdx-1)/2+2)) w l ls cs[cdx] lw mlw dt dts[cdx] t lbs[cdx]
set title 'M2 - P2' offset 0, -0.5
plot for [cdx=2:4:2] gx.'log' u ($0/1000):(column((cdx-2)/2+2)) w l ls cs[cdx] lw mlw dt dts[cdx] t lbs[cdx]
unset multiplot


set term wxt 5
set multiplot layout 1, 2
N = 180
set xrange [-0.5:N-0.5]
set yrange [N-0.5:-0.5]
set size square
set tics front
unset xlabel
unset ylabel
set xtics 45,90,135
set ytics 45,90,135
adx = 1
do for [cdx=1:1] {
  set arrow adx from graph 0, first cdx*90 to graph 1, first cdx*90 ls 8 front nohead
  adx = adx+1
  set arrow adx from first cdx*90, graph 0 to first cdx*90, graph 1 ls 8 front nohead
  adx = adx+1
}
load '../codes/reds.pal'
set logscale cb 10
set cbrange [0.005:1]
set title 'M1 - P1' offset 0, -0.5
do for [cdx=1:3:2] {
  set xtics add(lbs[cdx] cdx*45)
  set ytics add(lbs[cdx] cdx*45)
}
plot ux.'rc1.5_fq0.px' matrix u 1:2:3 w image notitle
set title 'M2 - P2' offset 0, -0.5
do for [cdx=2:4:2] {
  set xtics add(lbs[cdx] cdx*45-45)
  set ytics add(lbs[cdx] cdx*45-45)
}
plot vx.'rc1.5_fq0.px' matrix u 1:2:3 w image notitle
unset multiplot


