grep -i rugby_ball.log -e "norm" > norm
gnuplot -persist <<-EOFMarker
set logscale y 10
plot "norm" using 0:2
EOFMarker
rm norm
